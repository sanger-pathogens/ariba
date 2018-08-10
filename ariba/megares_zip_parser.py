import os
import sys
import csv
import zipfile
import pyfastaq
from ariba import common

class Error (Exception): pass


class MegaresZipParser:
    def __init__(self, zip_url, outprefix):
        self.zip_url = zip_url
        self.outprefix = outprefix
        self.zip_file = self.outprefix + '.downloaded.zip'


    @classmethod
    def _extract_files(cls, zip_file, outdir):
        original_files = {'annotations': None, 'fasta': None, 'header_mappings': None}

        try:
            os.mkdir(outdir)
        except:
            raise Error('Error making directory ' + outdir)

        zfile = zipfile.ZipFile(zip_file)
        for member in zfile.namelist():
            if '_annotations_' in member:
                original_files['annotations'] = member
            elif '_database_' in member and member.endswith('.fasta'):
                original_files['fasta'] = member
            elif '_header_mappings_' in member:
                original_files['header_mappings'] = member
            else:
                continue

            zfile.extract(member, path=outdir)

        if None in original_files.values():
            common.rmtree(outdir)
            raise Error('Error. Not all expected files found in downloaded megares zipfile. ' + str(original_files))

        return original_files


    @classmethod
    def _csv_to_dict(cls, infile, delimiter, expected_columns, key_column):
        data = {}
        non_key_columns = expected_columns - {key_column}

        with open(infile) as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            if set(expected_columns) != set(reader.fieldnames):
                raise Error('Unexpected header in annotations file. Expected columns: ' + ','.join(expected_columns) + ' but got: ' + ','.join(reader.fieldnames))

            for row in reader:
                data[row[key_column]] = {x: row[x] for x in non_key_columns}

        return data


    @classmethod
    def _load_annotations_file(cls, infile):
        return MegaresZipParser._csv_to_dict(infile, ',', {'header', 'class', 'mechanism', 'group'}, 'header')


    @classmethod
    def _load_header_mappings_file(cls, infile):
        return MegaresZipParser._csv_to_dict(infile, '\t', {'Source_Database', 'MEGARes_Header', 'Source_Headers(space_separated)'}, 'MEGARes_Header')


    @classmethod
    def _write_files(cls, outprefix, sequences, annotations, header_mappings):
        fasta = outprefix + '.fa'
        tsv = outprefix + '.tsv'
        fh_fasta = pyfastaq.utils.open_file_write(fasta)
        fh_tsv = pyfastaq.utils.open_file_write(tsv)

        for seq in sorted(sequences):
            final_column = []

            if seq in annotations:
                group = annotations[seq]['group']
                final_column.append('class:' + annotations[seq]['class'] + '; mechanism:' + annotations[seq]['mechanism'] + '; group:' + group)
            else:
                group = 'unknown'
                print('WARNING: sequence "', seq, '" has no record in annotations file', sep='', file=sys.stderr)

            if seq in header_mappings:
                final_column.append('Source_Database:' + header_mappings[seq]['Source_Database'] + '; Source_Headers:' + header_mappings[seq]['Source_Headers(space_separated)'])
            else:
                print('WARNING: sequence "', seq, '" has no record in header mappings file', sep='', file=sys.stderr)

            if len(final_column) > 0:
                print(group + '.' + seq, '1', '0', '.', '.', '; '.join(final_column), sep='\t', file=fh_tsv)
            else:
                print(group + '.' + seq, '1', '0', '.', '.', '.', sep='\t', file=fh_tsv)

            sequences[seq].id = group + '.' + sequences[seq].id
            print(sequences[seq], file=fh_fasta)

        fh_fasta.close()
        fh_tsv.close()


    def run(self):
        common.download_file(self.zip_url, self.zip_file, verbose=True)
        tmpdir = self.zip_file + '.tmp.extract'
        original_files = MegaresZipParser._extract_files(self.zip_file, tmpdir)
        annotation_data = MegaresZipParser._load_annotations_file(os.path.join(tmpdir, original_files['annotations']))
        header_data = MegaresZipParser._load_header_mappings_file(os.path.join(tmpdir, original_files['header_mappings']))
        sequences = {}
        pyfastaq.tasks.file_to_dict(os.path.join(tmpdir, original_files['fasta']), sequences)
        MegaresZipParser._write_files(self.outprefix, sequences, annotation_data, header_data)
        common.rmtree(tmpdir)
        os.unlink(self.zip_file)

