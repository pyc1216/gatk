package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Integration tests for {@link CNNScoreVariants}.
 * Created by sam on 1/8/18.
 */
public class CNNScoreVariantsIntegrationTest extends CommandLineProgramTest {
    private static final String architecture1D = largeFileTestDir + "VQSR/cnn_ref_model/1d_cnn_mix_train_full_bn.json";
    private static final String weights1D = largeFileTestDir + "VQSR/cnn_ref_model/1d_cnn_mix_train_full_bn.hd5";
    private static final String architecture2D = largeFileTestDir + "VQSR/cnn_read_model/small_2d.json";
    private static final String weights2D = largeFileTestDir + "VQSR/cnn_read_model/small_2d.hd5";
    private static final String inputVCF = largeFileTestDir + "VQSR/recalibrated_chr20_start.vcf";
    private static final String bigInputVCF = largeFileTestDir + "VQSR/g94982_20_1m_10m_python_2dcnn.vcf.gz";
    private static final String inputBAM = largeFileTestDir + "VQSR/g94982_contig_20_start_bamout.bam";
    private static final String inputIntervals = largeFileTestDir + "VQSR/contig20_conf_1m_10m.interval_list";
    private static final double EPSILON = 0.01;
    /**
     * Run the tool on a small test VCF.
     */
    @Test(groups = {"python"})
    public void testAllDefaultArgs() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                Arrays.asList(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf"));
        spec.executeTest("testInference", this);
    }

    /**
     * Run the tool on a small test VCF.
     */
    @Test(groups = {"python"})
    public void testAllDefaultArgsNew() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        final File tempVcf = createTempFile("tester", ".vcf");
        final File expectedVcf = new File(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf");
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath())
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        runCommandLine(argsBuilder);
        assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_1D_KEY);
    }

    @Test(groups = {"python"})
    public void testInferenceArchitecture() {
        final boolean newExpectations = false;
        final String expectedVCFName = largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("architecture", architecture1D)
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        if (newExpectations) {
            argsBuilder.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, expectedVCFName);
            runCommandLine(argsBuilder);
        } else {
            final File tempVcf = createTempFile("tester", ".vcf");
            final File expectedVcf = new File(expectedVCFName);
            argsBuilder.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath());
            runCommandLine(argsBuilder);
            assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_1D_KEY);
        }
    }

    @Test(groups = {"python"})
    public void testInferenceWeights() {
        final File tempVcf = createTempFile("tester", ".vcf");
        final File expectedVcf = new File(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("weights", weights1D)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath())
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        runCommandLine(argsBuilder);
        assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_1D_KEY);
    }

    @Test(groups = {"python"})
    public void testInferenceArchitectureAndWeights() {
        final File tempVcf = createTempFile("tester", ".vcf");
        final File expectedVcf = new File(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("weights", weights1D)
                .addArgument("architecture", architecture1D)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath())
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        runCommandLine(argsBuilder);
        assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_1D_KEY);
    }

    @Test(groups = {"python"})
    public void testInferenceWithIntervals() {
        final boolean newExpectations = false;
        final String expectedVCFName = largeFileTestDir + "VQSR/expected/cnn_1d_contig20_1m_10m_expected.vcf";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, bigInputVCF)
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument(StandardArgumentDefinitions.INTERVALS_LONG_NAME, inputIntervals)
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        if (newExpectations) {
            argsBuilder.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, expectedVCFName);
            runCommandLine(argsBuilder);
        } else {
            final File expectedVcf = new File(expectedVCFName);
            final File tempVcf = createTempFile("tester", ".vcf");
            argsBuilder.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath());
            runCommandLine(argsBuilder);
            assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_1D_KEY);
        }
    }

    @Test(groups = {"python"})
    public void testSmallBatchInference() {
        final File tempVcf = createTempFile("tester", ".vcf");
        final File expectedVcf = new File(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath())
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("inference-batch-size", "8")
                .addArgument("transfer-batch-size", "16")
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");
        runCommandLine(argsBuilder);
        assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_1D_KEY);
    }

    @Test(groups = {"python"})
    public void testOnContigEdge() {
        final String edgeVcf = toolsTestDir + "walkers/VQSR/variantNearContigEdge.vcf";
        final File tempVcf = createTempFile("tester", ".vcf");
        final File expectedVcf = new File(largeFileTestDir + "VQSR/expected/chrM.vcf");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, edgeVcf)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath())
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, hg19MiniReference)
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        runCommandLine(argsBuilder);
        assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_1D_KEY);
    }

    /**
     * Run the 2D Model on a small test VCF with the resource loaded weights and architecture.
     */
    @Test(groups = {"python"})
    public void testInference2dResourceModel() {
        // We reset the random number generator at the beginning of each test so that the random down-sampling of reads
        // by the reservoir down-sampler does not cause slightly different scores.
        Utils.resetRandomGenerator();
        TensorType tt = TensorType.read_tensor;
        final File tempVcf = createTempFile("tester", ".vcf");
        final File expectedVcf = new File(largeFileTestDir + "VQSR/expected/cnn_2d_chr20_subset_expected.vcf");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.INPUT_LONG_NAME, inputBAM)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath())
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("inference-batch-size", "2")
                .addArgument("transfer-batch-size", "2")
                .addArgument("tensor-type", tt.name())
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        runCommandLine(argsBuilder);
        assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_2D_KEY);
    }

    /**
     * Run the 2D Model on a small test VCF.
     */
    @Test(groups = {"python"})
    public void testInferenceArchitecture2d() {
        Utils.resetRandomGenerator();
        final boolean newExpectations = false;
        TensorType tt = TensorType.read_tensor;
        final String expectedVCFName = largeFileTestDir + "VQSR/expected/cnn_2d_chr20_subset_expected.vcf";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.INPUT_LONG_NAME, inputBAM)
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("architecture", architecture2D)
                .addArgument("tensor-type", tt.name())
                .addArgument("inference-batch-size", "8")
                .addArgument("transfer-batch-size", "8")
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        if (newExpectations) {
            argsBuilder.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, expectedVCFName);
            runCommandLine(argsBuilder);
        } else {
            final File tempVcf = createTempFile("tester", ".vcf");
            final File expectedVcf = new File(expectedVCFName);
            argsBuilder.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath());
            runCommandLine(argsBuilder);
            assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_2D_KEY);
        }
    }

    @Test(groups = {"python"})
    public void testInferenceWeights2d() {
        Utils.resetRandomGenerator();
        TensorType tt = TensorType.read_tensor;
        final File tempVcf = createTempFile("tester", ".vcf");
        final File expectedVcf = new File(largeFileTestDir + "VQSR/expected/cnn_2d_chr20_subset_expected.vcf");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.INPUT_LONG_NAME, inputBAM)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath())
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("weights", weights2D)
                .addArgument("inference-batch-size", "4")
                .addArgument("transfer-batch-size", "4")
                .addArgument("tensor-type", tt.name())
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        runCommandLine(argsBuilder);
        assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_2D_KEY);
    }

    @Test(groups = {"python"})
    public void testInferenceArchitectureAndWeights2d() {
        Utils.resetRandomGenerator();
        TensorType tt = TensorType.read_tensor;
        final File tempVcf = createTempFile("tester", ".vcf");
        final File expectedVcf = new File(largeFileTestDir + "VQSR/expected/cnn_2d_chr20_subset_expected.vcf");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.INPUT_LONG_NAME, inputBAM)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getPath())
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("weights", weights2D)
                .addArgument("architecture", architecture2D)
                .addArgument("inference-batch-size", "4")
                .addArgument("transfer-batch-size", "4")
                .addArgument("tensor-type", tt.name())
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        runCommandLine(argsBuilder);
        assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_2D_KEY);
    }

    private void assertInfoFieldsAreClose(File vcf1, File vcf2, String infoKey){
        Iterator<VariantContext> vi1 = VariantContextTestUtils.streamVcf(vcf1).collect(Collectors.toList()).iterator();
        Iterator<VariantContext> vi2 = VariantContextTestUtils.streamVcf(vcf2).collect(Collectors.toList()).iterator();
        while (vi1.hasNext() && vi2.hasNext()) {
            VariantContext vc1 = vi1.next();
            VariantContext vc2 = vi2.next();
            double v1Score = vc1.getAttributeAsDouble(infoKey, 0.0); // Different defaults trigger failures on missing scores
            double v2Score = vc2.getAttributeAsDouble(infoKey, EPSILON+1.0);
            double diff = Math.abs(v1Score-v2Score);
            Assert.assertTrue(diff < EPSILON);
        }
        Assert.assertTrue(!vi1.hasNext() && !vi2.hasNext());
    }



}
