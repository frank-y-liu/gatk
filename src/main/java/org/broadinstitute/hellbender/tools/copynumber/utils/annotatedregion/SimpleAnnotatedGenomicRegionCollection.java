package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hdf5.Utils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable.XsvLocatableTableCodec;
import org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable.XsvTableFeature;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.IntStream;

/**
 * Represents a collection of annotated regions.  The annotations do not need to be known ahead of time, if reading from a file.
 * Though {@link SimpleAnnotatedGenomicRegion::CONTIG_HEADER},
 *  {@link SimpleAnnotatedGenomicRegion::START_HEADER}, and
 *  {@link SimpleAnnotatedGenomicRegion::END_HEADER} are expected, by default, for defining the genomic region.
 */
public class SimpleAnnotatedGenomicRegionCollection {

    private SAMFileHeader samFileHeader;

    /** Does not include the locatable fields. */
    private List<String> annotations;
    private List<String> comments;
    private List<SimpleAnnotatedGenomicRegion> records;
    private String contigColumnName;
    private String startColumnName;
    private String endColumnName;

    private SimpleAnnotatedGenomicRegionCollection(final SAMFileHeader samFileHeader, final List<String> annotations, final List<String> comments,
                                                   final List<SimpleAnnotatedGenomicRegion> records, final String contigColumnName,
                                                   final String startColumnName, final String endColumnName) {
        this.comments = comments;
        this.samFileHeader = samFileHeader;
        this.annotations = annotations;
        this.records = records;

        /* The output contig column name */
        this.contigColumnName = contigColumnName;

        /* The output start column name */
        this.startColumnName = startColumnName;

        /* The output end column name */
        this.endColumnName = endColumnName;
    }

    //TODO: Docs
    public static SimpleAnnotatedGenomicRegionCollection create(final Path input, final Set<String> headersOfInterest) {
        ClassLoader classLoader = SimpleAnnotatedGenomicRegionCollection.class.getClassLoader();
        final File defaultConfigFile = new File(classLoader.getResource("org/broadinstitute/hellbender/tools/copynumber/utils/annotatedregion/annotated_region_default.config").getFile());
        return create(input, defaultConfigFile.toPath(), headersOfInterest);
    }

    // TODO: Docs
    public static SimpleAnnotatedGenomicRegionCollection create(final File input, final Set<String> headersOfInterest) {
        return create(input.toPath(), headersOfInterest);
    }

    //TODO: Docs
    public static SimpleAnnotatedGenomicRegionCollection create(final List<SimpleAnnotatedGenomicRegion> regions,
                                                                final SimpleAnnotatedGenomicRegionCollection collection) {

        return new SimpleAnnotatedGenomicRegionCollection(collection.getSamFileHeader(), collection.getAnnotations(), Collections.emptyList(), regions,
                collection.getContigColumnName(), collection.getStartColumnName(), collection.getEndColumnName());
    }
    // TODO:
    public static SimpleAnnotatedGenomicRegionCollection create(final List<SimpleAnnotatedGenomicRegion> regions,
                                                                                    final SAMFileHeader samFileHeader,
                                                                                    final List<String> annotations,
                                                                                    final String contigColumnName,
                                                                                    final String startColumnName,
                                                                                    final String endColumnName) {

        return new SimpleAnnotatedGenomicRegionCollection(samFileHeader, annotations, Collections.emptyList(), regions,
                contigColumnName, startColumnName, endColumnName);
    }

    //TODO: Docs
    public static SimpleAnnotatedGenomicRegionCollection create(final Path input, final Path inputConfigFile, final Set<String> headersOfInterest) {

        //TODO: Test for existence/viability of both files.

        final XsvLocatableTableCodec codec = new XsvLocatableTableCodec(inputConfigFile);
        final List<SimpleAnnotatedGenomicRegion> regions = new ArrayList<>();

        if (codec.canDecode(input.toString())) {
            try (final InputStream fileInputStream = Files.newInputStream(input)) {

                // Lots of scaffolding to do reading here:
                final AsciiLineReaderIterator lineReaderIterator = new AsciiLineReaderIterator(AsciiLineReader.from(fileInputStream));
                final List<String> header = codec.readActualHeader(lineReaderIterator);
                checkAllHeadersOfInterestPresent(headersOfInterest, header);

                final List<String> featureCols = codec.getHeaderWithoutLocationColumns();

                while (lineReaderIterator.hasNext()) {
                    final XsvTableFeature feature = codec.decode(lineReaderIterator.next());
                    if (feature == null) {
                        continue;
                    }

                    final List<String> featureValues = feature.getValuesWithoutLocationColumns();

                    final SortedMap<String, String> annotations = new TreeMap<>();
                    IntStream.range(0, featureCols.size()).boxed()
                            .filter(i -> (headersOfInterest == null) || headersOfInterest.contains(featureCols.get(i)))
                            .forEach(i -> annotations.put(featureCols.get(i), featureValues.get(i)));

                    regions.add(new SimpleAnnotatedGenomicRegion(
                            new SimpleInterval(feature.getContig(), feature.getStart(), feature.getEnd()),
                            annotations));
                }

                return new SimpleAnnotatedGenomicRegionCollection(codec.createSamFileHeader(), codec.getHeaderWithoutLocationColumns(),
                        codec.getComments(), regions, codec.getFinalContigColumn(), codec.getFinalStartColumn(), codec.getFinalEndColumn());

            }
            catch ( final FileNotFoundException ex ) {
                throw new GATKException("Error - could not find test file: " + input, ex);
            }
            catch ( final IOException ex ) {
                throw new GATKException("Error - IO problem with file " + input, ex);
            }
        }
        else {
            throw new UserException.BadInput("Could not parse xsv file.");
        }
    }

    private static void checkAllHeadersOfInterestPresent(final Set<String> headersOfInterest, final List<String> header) {
        if ((headersOfInterest != null) && !header.containsAll(headersOfInterest)) {
            final Set<String> unusedColumnsOfInterest = Sets.difference(new HashSet<>(headersOfInterest), new HashSet<>(header));
            if (unusedColumnsOfInterest.size() > 0) {
                final List<String> missingColumns = new ArrayList<>(unusedColumnsOfInterest);
                throw new UserException.BadInput("Some columns of interest specified by the user were not seen in any input files: " + StringUtils.join(missingColumns, ", "));
            }
        }
    }

    //TODO: Docs
    public void write(final File outputFile) {
        try (final FileWriter writer = new FileWriter(outputFile)) {
            writer.write(StringUtils.join(getComments(), "\n"));
            writer.write(getSamFileHeader().getSAMString());
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
        try (final RecordWriter recordWriter = new RecordWriter(new FileWriter(outputFile, true),
                contigColumnName, startColumnName, endColumnName, annotations)) {
            recordWriter.writeAllRecords(records);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    public SAMFileHeader getSamFileHeader() {
        return samFileHeader;
    }

    public List<String> getAnnotations() {
        return annotations;
    }

    public List<String> getComments() {
        return comments;
    }

    public List<SimpleAnnotatedGenomicRegion> getRecords() {
        return records;
    }

    public int size() {
        return getRecords().size();
    }

    public String getContigColumnName() {
        return contigColumnName;
    }

    public String getStartColumnName() {
        return startColumnName;
    }

    public String getEndColumnName() {
        return endColumnName;
    }

    private final class RecordWriter extends TableWriter<SimpleAnnotatedGenomicRegion> {
        private List<String> annotations;
        private String contigColumnName;
        private String startColumnName;
        private String endColumnName;

        private RecordWriter(final Writer writer, final String contigColumnName, final String startColumnName,
                             final String endColumnName, final List<String> annotationsToWrite) throws IOException {
            super(writer, new TableColumnCollection(annotationsToWrite));
            this.annotations = annotationsToWrite;
            this.contigColumnName = contigColumnName;
            this.startColumnName = startColumnName;
            this.endColumnName = endColumnName;
        }

        @Override
        protected void composeLine(final SimpleAnnotatedGenomicRegion record, final DataLine dataLine) {
            Utils.nonNull(record);
            Utils.nonNull(dataLine);

            // First handle the locatable columns.
            dataLine.set(contigColumnName, record.getContig());
            dataLine.set(startColumnName, record.getStart());
            dataLine.set(endColumnName, record.getEnd());

            // Then the annotations.
            annotations.forEach(a -> dataLine.set(a, record.getAnnotationValue(a)));
        }
    }
}
