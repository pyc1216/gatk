package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.AnalyzeSaturationMutagenesis.Interval;
import org.broadinstitute.hellbender.tools.AnalyzeSaturationMutagenesis.IntervalCounter;
import org.broadinstitute.hellbender.tools.AnalyzeSaturationMutagenesis.SNV;
import org.broadinstitute.hellbender.tools.AnalyzeSaturationMutagenesis.SNVCollectionCount;
import org.broadinstitute.hellbender.tools.AnalyzeSaturationMutagenesis.CodonTracker;
import org.broadinstitute.hellbender.tools.AnalyzeSaturationMutagenesis.CodonVariation;
import org.broadinstitute.hellbender.tools.AnalyzeSaturationMutagenesis.CodonVariationType;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class AnalyzeSaturationMutagenesisUnitTest extends GATKBaseTest {
    @Test
    public void testInterval() {
        final Interval interval = new Interval(1, 3);
        Assert.assertEquals(interval.size(), 2);
    }

    @Test
    public void testSNV() {
        final SNV snv1 = new SNV(0, (byte)'A', (byte)'C', (byte)30);
        final SNV snv2 = new SNV(0, (byte)'A', (byte)'C', (byte)20);
        Assert.assertEquals(snv1.hashCode(), snv2.hashCode());
        Assert.assertEquals(snv1, snv2);
        Assert.assertEquals(snv1.compareTo(snv2), 0);
        final SNV snv3 = new SNV(1, (byte)'A', (byte)'C', (byte)30);
        Assert.assertNotEquals(snv1.hashCode(), snv3.hashCode());
        Assert.assertNotEquals(snv1, snv3);
        Assert.assertTrue(snv1.compareTo(snv3) < 0);
        Assert.assertTrue(snv3.compareTo(snv2) > 0);
        final SNV snv4 = new SNV(0, (byte)'G', (byte)'C', (byte)30);
        Assert.assertNotEquals(snv1.hashCode(), snv4.hashCode());
        Assert.assertNotEquals(snv1, snv4);
        Assert.assertTrue(snv1.compareTo(snv4) < 0);
        Assert.assertTrue(snv4.compareTo(snv1) > 0);
        final SNV snv5 = new SNV(0, (byte)'A', (byte)'G', (byte)30);
        Assert.assertNotEquals(snv1.hashCode(), snv5.hashCode());
        Assert.assertNotEquals(snv1, snv5);
        Assert.assertTrue(snv1.compareTo(snv5) < 0);
        Assert.assertTrue(snv5.compareTo(snv1) > 0);
    }

    @Test
    public void testIntervalCounter() {
        final IntervalCounter intervalCounter = new IntervalCounter(10);
        intervalCounter.addCount(0, 10);
        intervalCounter.addCount(1, 9);
        intervalCounter.addCount(2, 8);
        intervalCounter.addCount(3, 7);
        intervalCounter.addCount(4, 6);
        Assert.assertEquals(intervalCounter.countSpanners(0, 10), 1);
        Assert.assertEquals(intervalCounter.countSpanners(5, 5), 5);
        Assert.assertEquals(intervalCounter.countSpanners(2, 5), 3);
        Assert.assertEquals(intervalCounter.countSpanners(5, 8), 3);
        intervalCounter.addCount(0, 10);
        Assert.assertEquals(intervalCounter.countSpanners(0, 10), 2);
    }

    @Test
    public void testSNVCollection() {
        final List<SNV> snvList = new ArrayList<>(Arrays.asList(
                new SNV(0, (byte)'A', (byte)'C', (byte)30),
                new SNV(1, (byte)'A', (byte)'C', (byte)20)));
        final SNVCollectionCount cc1 = new SNVCollectionCount(snvList, 10);

        // equality, key, compare, and hash should be independent of count and coverage
        final List<SNV> snvList2 = Arrays.asList(
                new SNV(0, (byte)'A', (byte)'C', (byte)30),
                new SNV(1, (byte)'A', (byte)'C', (byte)20));
        final SNVCollectionCount cc2 = new SNVCollectionCount(snvList2, 20);
        Assert.assertEquals(cc1.hashCode(), cc2.hashCode());
        Assert.assertEquals(cc1, cc2);
        Assert.assertEquals(cc1.compareTo(cc2), 0);
        cc2.bumpCount(30);
        Assert.assertEquals(cc1.hashCode(), cc2.hashCode());
        Assert.assertEquals(cc1, cc2);
        Assert.assertEquals(cc1.compareTo(cc2), 0);

        Assert.assertEquals(cc2.getCount(), 2);
        Assert.assertEquals(cc2.getMeanRefCoverage(), 25., .0000001);

        // changing the list shouldn't change the hash or the key
        final int cc1Hash = cc1.hashCode();
        final List<SNV> key1 = cc1.getSNVs();
        snvList.add(new SNV(2, (byte)'A', (byte)'C', (byte)10));
        Assert.assertEquals(cc1.hashCode(), cc1Hash);
        Assert.assertEquals(cc1.getSNVs(), key1);

        // different lists should mean unequal to each other, unequal hashes, and non-zero compare
        final SNVCollectionCount cc3 = new SNVCollectionCount(snvList, 20);
        Assert.assertNotEquals(cc1.hashCode(), cc3.hashCode());
        Assert.assertNotEquals(cc1, cc3);
        Assert.assertTrue(cc1.compareTo(cc3) < 0);
        Assert.assertTrue(cc3.compareTo(cc1) > 0);
    }

    @Test
    public void testCodonVariation() {
        Assert.assertTrue(new CodonVariation(0,0, CodonVariationType.DELETION).isDeletion());
        Assert.assertTrue(new CodonVariation(0,0, CodonVariationType.FRAMESHIFT).isFrameshift());
        Assert.assertTrue(new CodonVariation(0,0, CodonVariationType.INSERTION).isInsertion());
        Assert.assertTrue(new CodonVariation(0,0, CodonVariationType.MODIFICATION).isModification());
        Assert.assertFalse(new CodonVariation(0,0, CodonVariationType.DELETION).isModification());
    }

    final static byte CALL_A = (byte)'A';
    final static byte CALL_C = (byte)'C';
    final static byte CALL_G = (byte)'G';
    final static byte CALL_T = (byte)'T';
    final static byte NO_CALL = (byte)'-';
    final static byte QUAL_30 = (byte)30;

    final static byte[] refSeq = "ACATGCGTCTAGTACGT".getBytes();
    final static String orfCoords = "3-6,8-12";
    final static CodonTracker codonTracker = new CodonTracker(orfCoords, refSeq);

    @Test
    public void testEncodingModifications() {
        // no SNVs implies no codon variations
        Assert.assertTrue(codonTracker.encodeSNVsAsCodons(Collections.emptyList()).isEmpty());

        // changes outside the ORF shouldn't produce codon variations
        Assert.assertTrue(codonTracker
                .encodeSNVsAsCodons(Collections.singletonList(new SNV(1, CALL_C, CALL_G, QUAL_30)))
                .isEmpty());
        Assert.assertTrue(codonTracker
                .encodeSNVsAsCodons(Collections.singletonList(new SNV(6, CALL_G, CALL_A, QUAL_30)))
                .isEmpty());
        Assert.assertTrue(codonTracker
                .encodeSNVsAsCodons(Collections.singletonList(new SNV(12, CALL_T, CALL_C, QUAL_30)))
                .isEmpty());

        // changes to a single codon should produce a single-codon variations
        Assert.assertEquals(
                codonTracker.encodeSNVsAsCodons(Collections.singletonList(new SNV(2, CALL_A, CALL_C, QUAL_30))),
                Collections.singletonList(new CodonVariation(0, 30, CodonVariationType.MODIFICATION)));
        Assert.assertEquals(
                codonTracker.encodeSNVsAsCodons(Collections.singletonList(new SNV(5, CALL_C, CALL_G, QUAL_30))),
                Collections.singletonList(new CodonVariation(1, 45, CodonVariationType.MODIFICATION)));
        Assert.assertEquals(
                codonTracker.encodeSNVsAsCodons(Collections.singletonList(new SNV(7, CALL_T, CALL_G, QUAL_30))),
                Collections.singletonList(new CodonVariation(1, 25, CodonVariationType.MODIFICATION)));
        Assert.assertEquals(
                codonTracker.encodeSNVsAsCodons(Collections.singletonList(new SNV(8, CALL_C, CALL_A, QUAL_30))),
                Collections.singletonList(new CodonVariation(1, 28, CodonVariationType.MODIFICATION)));
        Assert.assertEquals(
                codonTracker.encodeSNVsAsCodons(Collections.singletonList(new SNV(11, CALL_G, CALL_A, QUAL_30))),
                Collections.singletonList(new CodonVariation(2, 48, CodonVariationType.MODIFICATION)));
        Assert.assertEquals(
                codonTracker.encodeSNVsAsCodons(Arrays.asList(
                        new SNV(5, CALL_C, CALL_G, QUAL_30),
                        new SNV(7, CALL_T, CALL_G, QUAL_30),
                        new SNV(8, CALL_C, CALL_A, QUAL_30))),
                Collections.singletonList(new CodonVariation(1, 40, CodonVariationType.MODIFICATION)));

        // even if the change produces a nonsense codon
        Assert.assertEquals(
                codonTracker.encodeSNVsAsCodons(Arrays.asList(
                        new SNV(5, CALL_C, CALL_T, QUAL_30),
                        new SNV(7, CALL_T, CALL_A, QUAL_30),
                        new SNV(8, CALL_C, CALL_A, QUAL_30))),
                Collections.singletonList(new CodonVariation(1, 48, CodonVariationType.MODIFICATION)));
    }

    @Test
    public void testEncodingDeletions() {
        Assert.assertEquals(
                codonTracker.encodeSNVsAsCodons(Collections.singletonList(new SNV(2, CALL_A, NO_CALL, QUAL_30))),
                Arrays.asList(
                        new CodonVariation(0, -1, CodonVariationType.FRAMESHIFT),
                        new CodonVariation(0, 57, CodonVariationType.MODIFICATION),
                        new CodonVariation(1, 55, CodonVariationType.MODIFICATION),
                        new CodonVariation(2, 11, CodonVariationType.MODIFICATION)));
        Assert.assertEquals(
                codonTracker.encodeSNVsAsCodons(Arrays.asList(
                        new SNV(10, CALL_A, NO_CALL, QUAL_30),
                        new SNV(12, CALL_T, NO_CALL, QUAL_30))),
                Arrays.asList(
                        new CodonVariation(2, -1, CodonVariationType.FRAMESHIFT),
                        new CodonVariation(2, 56, CodonVariationType.MODIFICATION)));
        Assert.assertEquals(
                codonTracker.encodeSNVsAsCodons(Arrays.asList(
                        new SNV(5, CALL_C, NO_CALL, QUAL_30),
                        new SNV(7, CALL_T, NO_CALL, QUAL_30),
                        new SNV(8, CALL_C, NO_CALL, QUAL_30))),
                Collections.singletonList(new CodonVariation(1, -1, CodonVariationType.DELETION)));
        Assert.assertEquals(
                codonTracker.encodeSNVsAsCodons(Arrays.asList(
                        new SNV(3, CALL_T, NO_CALL, QUAL_30),
                        new SNV(5, CALL_C, NO_CALL, QUAL_30),
                        new SNV(8, CALL_C, NO_CALL, QUAL_30))),
                Arrays.asList(
                        new CodonVariation(0, 11, CodonVariationType.MODIFICATION),
                        new CodonVariation(1, -1, CodonVariationType.DELETION)));

    }

    @Test
    public void testEncodingInsertions() {
        Assert.assertEquals(
                codonTracker.encodeSNVsAsCodons(Collections.singletonList(
                        new SNV(5, NO_CALL, CALL_T, QUAL_30))),
                Arrays.asList(
                        new CodonVariation(1, -1, CodonVariationType.FRAMESHIFT),
                        new CodonVariation(1, 55, CodonVariationType.MODIFICATION),
                        new CodonVariation(2, 28, CodonVariationType.MODIFICATION)));
        Assert.assertEquals(
                codonTracker.encodeSNVsAsCodons(Arrays.asList(
                        new SNV(2, NO_CALL, CALL_T, QUAL_30),
                        new SNV(2, NO_CALL, CALL_T, QUAL_30),
                        new SNV(2, NO_CALL, CALL_T, QUAL_30))),
                Collections.singletonList(
                        new CodonVariation(0, 63, CodonVariationType.INSERTION)));
    }

    @Test
    public void testFrameRecoveringIndels() {
        Assert.assertEquals(
                codonTracker.encodeSNVsAsCodons(Arrays.asList(
                        new SNV(3, CALL_T, NO_CALL, QUAL_30),
                        new SNV(7, NO_CALL, CALL_A, QUAL_30))),
                Arrays.asList(
                        new CodonVariation(0, 9, CodonVariationType.MODIFICATION),
                        new CodonVariation(1, 13, CodonVariationType.MODIFICATION)));
        Assert.assertEquals(
                codonTracker.encodeSNVsAsCodons(Arrays.asList(
                        new SNV(3, NO_CALL, CALL_T, QUAL_30),
                        new SNV(7, CALL_T, NO_CALL, QUAL_30))),
                Arrays.asList(
                        new CodonVariation(0, 15, CodonVariationType.MODIFICATION),
                        new CodonVariation(1, 37, CodonVariationType.MODIFICATION)));
    }
}
