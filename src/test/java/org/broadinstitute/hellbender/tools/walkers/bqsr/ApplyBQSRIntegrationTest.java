package org.broadinstitute.hellbender.tools.walkers.bqsr;

import htsjdk.samtools.SamReaderFactory;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;

public final class ApplyBQSRIntegrationTest extends CommandLineProgramTest {
    private static class ABQSRTest {
        final String bam;
        final String reference;
        final String outputExtension;
        final String args[];
        final String expectedFile;

        private ABQSRTest(String bam, String reference, String outputExtension, String args[], String expectedFile) {
            this.bam= bam;
            this.reference = reference;
            this.outputExtension = outputExtension;
            this.args = args;
            this.expectedFile = expectedFile;
        }

        @Override
        public String toString() {
            return String.format("ApplyBQSR(args='%s')", args == null ? "" : StringUtils.join(args));
        }
    }

    @Override
    public String getTestedClassName() {
        return ApplyBQSR.class.getSimpleName();
    }

    final String resourceDir = getTestDataDir() + "/" + "BQSR" + "/";
    final String hg18Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";
    final String hiSeqBam = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate.bam";
    final String hiSeqCram = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate.cram";
    final String hiSeqBamAligned = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate_allaligned.bam";
    final String hiSeqCramAligned = resourceDir + "HiSeq.1mb.1RG.2k_lines.alternate_allaligned.cram";

    @DataProvider(name = "ApplyBQSRTest")
    public Object[][] createABQSRTestData() {
        List<Object[]> tests = new ArrayList<>();

        //Note: these outputs were created using GATK3
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", null, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", new String[] {"-OQ"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.OQ.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", new String[] {"-qq", "-1"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.qq-1.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", new String[] {"-qq", "6"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.qq6.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", new String[] {"-SQQ", "10", "-SQQ", "20", "-SQQ", "30"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.SQQ102030.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBam, null, ".bam", new String[] {"-SQQ", "10", "-SQQ", "20", "-SQQ", "30", "-RDQ"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.SQQ102030RDQ.bam")});

        tests.add(new Object[]{new ABQSRTest(hiSeqBamAligned, null, ".bam", null, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBamAligned, null, ".bam", new String[] {"-OQ"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.OQ.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBamAligned, null, ".bam", new String[] {"-qq", "-1"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.qq-1.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBamAligned, null, ".bam", new String[] {"-qq", "6"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.qq6.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBamAligned, null, ".bam", new String[] {"-SQQ", "10", "-SQQ", "20", "-SQQ", "30"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.SQQ102030.bam")});
        tests.add(new Object[]{new ABQSRTest(hiSeqBamAligned, null, ".bam", new String[] {"-SQQ", "10", "-SQQ", "20", "-SQQ", "30", "-RDQ"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.SQQ102030RDQ.bam")});

        //CRAM - input and output crams generated by direct conversion of the corresponding BAM test files with samtools 1.3
        tests.add(new Object[]{new ABQSRTest(hiSeqCram, hg18Reference, ".cram", new String[] {"--disableSequenceDictionaryValidation", "true"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate.recalibrated.DIQ.cram")});
        tests.add(new Object[]{new ABQSRTest(hiSeqCramAligned, hg18Reference, ".cram", new String[] {"-qq", "6", "--disableSequenceDictionaryValidation", "true"}, resourceDir + "expected.HiSeq.1mb.1RG.2k_lines.alternate_allaligned.recalibrated.DIQ.qq6.cram")});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ApplyBQSRTest")
    public void testApplyBQSR(ABQSRTest params) throws IOException {
        File outFile = BaseTest.createTempFile("applyBQSRTest", params.outputExtension);
        final ArrayList<String> args = new ArrayList<>();
        File refFile = null;

        args.add("-I");
        args.add(new File(params.bam).getAbsolutePath());
        args.add("--bqsr_recal_file");
        args.add(new File(resourceDir + "HiSeq.20mb.1RG.table.gz").getAbsolutePath());
        args.add("-O"); args.add(outFile.getAbsolutePath());
        if (params.reference != null) {
            refFile = new File(params.reference);
            args.add("-R"); args.add(refFile.getAbsolutePath());
        }
        if (params.args != null) {
            Stream.of(params.args).forEach(arg -> args.add(arg));
        }

        runCommandLine(args);

        SamAssertionUtils.assertSamsEqual(outFile, new File(params.expectedFile), refFile);
    }

    @Test
    public void testMissingReadGroup() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -I " + hiSeqBamAligned +
                        " --bqsr_recal_file " + resourceDir + "HiSeq.20mb.1RG.table.missingRG.gz" +
                        " -O /dev/null", 0,
                IllegalStateException.class);
        spec.executeTest("testMissingReadGroup", this);
    }

    @Test
    public void testemptyBqsrRecalFile() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -I " + hiSeqBamAligned +
                        " --bqsr_recal_file " + createTempFile("emptyBqsrRecal", "").toString() +
                        " -O /dev/null", 0,
                UserException.class);
        spec.executeTest("testemptyBqsrRecalFile", this);
    }

    @Test
    public void testPRNoFailWithHighMaxCycle() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -I " + hiSeqBamAligned +
                        " --bqsr_recal_file " + resourceDir + "HiSeq.1mb.1RG.highMaxCycle.table.gz" +
                        " -O /dev/null",
                Arrays.<String>asList());
        spec.executeTest("testPRNoFailWithHighMaxCycle", this);      //this just checks that the tool does not blow up
    }


    @Test
    public void testHelp() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -I " + hiSeqBamAligned +
                        " --help --bqsr_recal_file " + resourceDir + "HiSeq.1mb.1RG.highMaxCycle.table.gz" +
                        " -O /dev/null",
                Arrays.<String>asList());
        spec.executeTest("testHelp", this);      //this just checks that the tool does not blow up
    }

    @Test
    public void testPRFailWithLowMaxCycle() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -I " + hiSeqBamAligned +
                        " --bqsr_recal_file " + resourceDir + "HiSeq.1mb.1RG.lowMaxCycle.table.gz" +
                        " -O /dev/null",
                0,
                UserException.class);
        spec.executeTest("testPRFailWithLowMaxCycle", this);
    }

    @Test
    public void testPRWithConflictingArguments_qqAndSQQ() throws IOException {
        // -qq and -SQQ shouldn't be able to be run in the same command
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -I " + hiSeqBam +
                        " --bqsr_recal_file " + resourceDir + "HiSeq.20mb.1RG.table.gz" +
                        " -SQQ 9 -qq 4 " +
                        " -O /dev/null",
                0,
                CommandLineException.class);
        spec.executeTest("testPRWithConflictingArguments_qqAndSQQ", this);
    }

    @Test
    public void testPRWithConflictingArguments_qqAndRDQ() throws IOException {
        // -qq and -SQQ shouldn't be able to be run in the same command
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -I " + hiSeqBam +
                        " --bqsr_recal_file " + resourceDir + "HiSeq.20mb.1RG.table.gz" +
                        " -RDQ -qq 4 " +
                        " -O /dev/null",
                0,
                CommandLineException.class);
        spec.executeTest("testPRWithConflictingArguments_qqAndSQQ", this);
    }

    @Test
    public void testOverfiltering() throws IOException {
        final File zeroRefBasesReadBam = new File(resourceDir, "NA12878.oq.read_consumes_zero_ref_bases.bam");
        final File outFile = BaseTest.createTempFile("testReadThatConsumesNoReferenceBases", ".bam");
        final String[] args = new String[] {
                "--input", zeroRefBasesReadBam.getAbsolutePath(),
                "--bqsr_recal_file", resourceDir + "NA12878.oq.gatk4.recal.gz",
                "--useOriginalQualities",
                "--output", outFile.getAbsolutePath()
        };
        runCommandLine(args);
        //The expected output is actually the same as inputs for this read
        SamAssertionUtils.assertSamsEqual(outFile, zeroRefBasesReadBam);
    }

    @Test
    public void testAddingPG() throws IOException {
        final File inFile = new File(resourceDir, "NA12878.oq.read_consumes_zero_ref_bases.bam");
        final File outFile = BaseTest.createTempFile("testAddingPG", ".bam");
        final String[] args = new String[] {
                "--input", inFile.getAbsolutePath(),
                "--bqsr_recal_file", resourceDir + "NA12878.oq.gatk4.recal.gz",
                "--useOriginalQualities",
                "--addOutputSAMProgramRecord",
                "--output", outFile.getAbsolutePath()
        };
        runCommandLine(args);

        //The expected output is actually the same as inputs for this read (this ignores the PGs header)
        SamAssertionUtils.assertSamsEqual(outFile, inFile);

        //input has no GATK ApplyBQSR in headers
        Assert.assertNull(SamReaderFactory.makeDefault().open(inFile).getFileHeader().getProgramRecord("GATK ApplyBQSR"));

        //output has a GATK ApplyBQSR in headers
        Assert.assertNotNull(SamReaderFactory.makeDefault().open(outFile).getFileHeader().getProgramRecord("GATK ApplyBQSR"));
    }
}
