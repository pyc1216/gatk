package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import static java.lang.Math.max;
import static java.lang.Math.min;

/**
 * Strand bias estimated by the Symmetric Odds Ratio test
 *
 * <p>Strand bias is a type of sequencing bias in which one DNA strand is favored over the other, which can result in
 * incorrect evaluation of the amount of evidence observed for one allele vs. the other. The StrandOddsRatio annotation
 * is one of several methods that aims to evaluate whether there is strand bias in the data. It is an updated form of
 * the <a href="https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_annotator_FisherStrand.php">Fisher Strand Test</a>
 * that is better at taking into account large amounts of data in high coverage situations. It is used to determine if
 * there is strand bias between forward and reverse strands for the reference or alternate allele(s).</p>
 *
 * <h3>Statistical notes</h3>
 *
 * <p>The following 2x2 contingency table gives the notation for allele support and strand orientation.</p>
 *
 * <table>
 * <tr><th>&nbsp;</th><th>+ strand&nbsp;&nbsp;&nbsp;</th><th>- strand&nbsp;&nbsp;&nbsp;</th></tr>
 * <tr><th>REF&nbsp;&nbsp;&nbsp;</th><td>X[0][0]</td><td>X[0][1]</td></tr>
 * <tr><th>ALT&nbsp;&nbsp;&nbsp;</th><td>X[1][0]</td><td>X[1][1]</td></tr>
 * </table>
 *
 * <p>We can then represent the Odds Ratios with the equation:</p>
 *
 * <img src="http://latex.codecogs.com/svg.latex?$$ R = \frac{X[0][0] * X[1][1]}{X[0][1] * X[1][0]} $$" border="0"/>
 *
 * <p>and its inverse:</p>
 *
 * <img src="http://latex.codecogs.com/svg.latex?$$ \frac{1}{R} = \frac{X[0][1] * X[1][0]}{X[0][0] * X[1][1]} $$" border="0"/>
 *
 * <p>The sum R + 1/R is used to detect a difference in strand bias for REF and for ALT. The sum makes it symmetric.
 * A high value is indicative of large difference where one entry is very small compared to the others. A scale factor
 * of refRatio/altRatio where</p>
 *
 * <img src="http://latex.codecogs.com/svg.latex?$$ refRatio = \frac{min(X[0][0], X[0][1])}{max(X[0][0], X[0][1])} $$" border="0"/>
 *
 * <p>and </p>
 *
 * <img src="http://latex.codecogs.com/svg.latex?$$ altRatio = \frac{min(X[1][0], X[1][1])}{max(X[1][0], X[1][1])} $$" border="0"/>
 *
 * <p>ensures that the annotation value is large only. The final SOR annotation is given in natural log space.</p>
 *
 * <p>See the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4732">method document on statistical tests</a>
 * for a more detailed explanation of this statistical test, and see
 * <a href="https://gatkforums.broadinstitute.org/gatk/discussion/5533/strandoddsratio-computation">this post</a> for an
 * example calculation.</p>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_annotator_allelespecific_AS_StrandOddsRatio.php">AS_StrandOddsRatio</a></b>
 *     allele-specific strand bias estimated by the symmetric odds ratio test.</li>
 *     <li><b><a href="https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_annotator_StrandBiasBySample.php">StrandBiasBySample</a></b>
 *     outputs counts of read depth per allele for each strand orientation.</li>
 *     <li><b><a href="https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_annotator_FisherStrand.php">FisherStrand</a></b>
 *     uses Fisher's Exact Test to evaluate strand bias.</li>
 * </ul>
 *
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Strand bias estimated by the symmetric odds ratio test (SOR)")
public final class StrandOddsRatio extends StrandBiasTest implements StandardAnnotation {

    private static final double PSEUDOCOUNT = 1.0;
    private static final int MIN_COUNT = 0;

    @Override
    protected Map<String, Object> calculateAnnotationFromGTfield(final GenotypesContext genotypes){
        final int[][] tableFromPerSampleAnnotations = getTableFromSamples(genotypes, MIN_COUNT);
        return tableFromPerSampleAnnotations != null ? annotationForOneTable(calculateSOR(tableFromPerSampleAnnotations)) : null;
    }

    @Override
    protected Map<String, Object> calculateAnnotationFromStratifiedContexts(Map<String, List<PileupElement>> stratifiedContexts,
                                                                            final VariantContext vc){
        final int[][] tableNoFiltering = getPileupContingencyTable(stratifiedContexts, vc.getReference(), vc.getAlternateAlleles(), -1, MIN_COUNT);
        final double ratio = calculateSOR(tableNoFiltering);
        return annotationForOneTable(ratio);
    }


    @Override
    protected Map<String, Object> calculateAnnotationFromLikelihoods(final ReadLikelihoods<Allele> likelihoods, final VariantContext vc){
        final int[][] table = getContingencyTable(likelihoods, vc, MIN_COUNT);
        return annotationForOneTable(calculateSOR(table));
    }

    /**
     * Computes the SOR value of a table after augmentation. Based on the symmetric odds ratio but modified to take on
     * low values when the reference +/- read count ratio is skewed but the alt count ratio is not.  Natural log is taken
     * to keep values within roughly the same range as other annotations.
     *
     * Adding pseudocounts avoids division by zero.
     *
     * @param table The table before adding pseudocounts
     * @return the SOR annotation value
     */
    public static double calculateSOR(final int[][] table) {
        final double t00 = table[0][0] + PSEUDOCOUNT;
        final double t01 = table[0][1] + PSEUDOCOUNT;
        final double t11 = table[1][1] + PSEUDOCOUNT;
        final double t10 = table[1][0] + PSEUDOCOUNT;

        final double ratio = (t00 / t01) * (t11 / t10) + (t01 / t00) * (t10 / t11);

        final double refRatio = min(t00, t01)/ max(t00, t01);
        final double altRatio = min(t10, t11)/ max(t10, t11);

        return Math.log(ratio) + Math.log(refRatio) - Math.log(altRatio);
    }

    /**
     * Returns an annotation result given a sor
     *
     * @param sor the symmetric odds ratio of the contingency table
     * @return a hash map from SOR
     */
    @VisibleForTesting
    Map<String, Object> annotationForOneTable(final double sor) {
        return Collections.singletonMap(getKeyNames().get(0), formattedValue(sor));
    }

    public static String formattedValue(double sor) {
        return String.format("%.3f", sor);
    }


    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GATKVCFConstants.STRAND_ODDS_RATIO_KEY);
    }
}
