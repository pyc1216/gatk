package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.testng.Assert.*;

/**
 * Created by gauthier on 6/4/18.
 */
public class VariantsToTableUnitTest {
    final Allele refAllele = Allele.create("A", true);
    final Allele altAllele1 = Allele.create("T");
    final Allele altAllele2 = Allele.create("C");

    @Test
    public void testExtractFields() throws Exception {
        final VariantsToTable VTTtest = VariantsToTable.getInstance();
        VTTtest.fieldsToTake = Arrays.asList("AC");
        VTTtest.splitMultiAllelic = true;

        final List<Allele> alleles = Arrays.asList(refAllele, altAllele1, altAllele2);
        final VariantContext vc = new VariantContextBuilder("test", "1", 10000, 10000, alleles).attribute(VCFConstants.ALLELE_COUNT_KEY, "1,3").make();
        final List<List<String>> fields = VariantsToTable.getInstance().extractFields(vc);
    }

    @Test
    public void testExtractRawAnnotations() {
        final Genotype genotype = new GenotypeBuilder("sample2", Arrays.asList(refAllele, altAllele1, altAllele2))
                .AD(new int[]{2,80,9}).make();

        final int MQsquaredRef = 400;
        final int MQsquaredAlt1 = 285;
        final int MQsquaredAlt2 = 385;

        final VariantContext vc = new VariantContextBuilder(new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele1, altAllele2))
                .chr("1").start(15L).stop(15L)
                .attribute(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY, "400.00|285.00|385.00")
                .genotypes(genotype)
                .make();

    }

}