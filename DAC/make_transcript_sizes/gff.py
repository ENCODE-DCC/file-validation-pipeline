# Copyright (c) 2007
# Colin Dewey (University of Wisconsin-Madison)
# cdewey@biostat.wisc.edu
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# References for GFF formats
# version   reference
# 1 & 2     http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml
# 2.1       http://genes.cs.wustl.edu/GTF21.html
# 2.2       http://genes.cs.wustl.edu/GTF22.html
# 2.5       ???
# 3         http://song.sourceforge.net/gff3.shtml

import sys
import re
from urllib import quote as url_quote, unquote as url_unquote
from collections import defaultdict

__all__ = ["Record", "Reader", "Writer", "FormatError",
           "Metadatum", "SequenceRegion"]

class Record:
    """A record from a GFF file.

    Fields:
       seqid
       source
       type
       start
       end
       score
       strand
       phase
       attributes
    """

    def __init__(self, seqid, source, type, start, end,
                 score=None, strand=None, phase=None, attributes=None):
        self.seqid = seqid
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase

        if attributes:
            self.attributes = attributes
        else:
            self.attributes = {}

    def copy(self):
        """Returns a copy of this GFF record"""

        # Make new attributes dict with copied value lists
        attributes_copy = dict([(k, v[:]) for k, v in self.attributes.items()])

        return Record(self.seqid,
                      self.source,
                      self.type,
                      self.start,
                      self.end,
                      score=self.score,
                      strand=self.strand,
                      phase=self.phase,
                      attributes=attributes_copy)

    def is_valid(self):
        """Returns True if this record passes basic GFF record requirements."""

        # Check that the first five fields are defined and non-empty/zero
        if not (self.seqid and self.source and self.type and
                self.start and self.end):
            return False

        # Check that [start, end] is a valid genomic interval
        if not (is_integer(self.start) and is_integer(self.end) and
                  self.start > 0 and self.start <= self.end):
            return False
        
        # Check that score is a valid floating point number
        if self.score is not None:
            try:
                float(self.score)
            except ValueError:
                return False
        
        if self.strand not in (None, '+', '-'):
            return False
        if self.phase not in (None, 0, 1, 2):
            return False
        
        # Check that CDS records have phase defined
        if self.type == "CDS" and self.phase is None:
            return False

        return True

    def __repr__(self):
        return "Record(%s, %s, %s, %s, %s, %s, %s, %s, %s)" % \
            tuple(map(repr, (self.seqid, self.source, self.type, 
                             self.start, self.end,
                             self.score, self.strand, self.phase, 
                             self.attributes)))

class Metadatum:
    def __init__(self, name, value=None):
        self.name = name
        self.value = value

class SequenceRegion(Metadatum):
    def __init__(self, seqid, start, end):
        Metadatum.__init__(self, 
                           "sequence-region", 
                           "%s %d %d" % (seqid, start, end))
        self.seqid = seqid
        self.start = start
        self.end = end

class FormatError(Exception):
    """Invalid format for GFF file"""
    pass

class Reader:
    """Reads a GFF formatted file"""

    def __init__(self, stream, version="2"):
        self._stream = stream
        self._default_version = version

        # Directives
        self._version = None
        self._references_resolved = True

        # Metadata
        self._metadata = []
        self._comments = []
        self._sequence_regions = []
        self._fasta_string = ""

        # Create record parser table
        self._record_parsers = {"1": self._parse_record_v1,
                                "2": self._parse_record_v2,
                                "2.1": self._parse_record_v2,
                                "2.2": self._parse_record_v2,
                                "2.5": self._parse_record_v2,
                                "3": self._parse_record_v3}

        #  Set default record parser
        if version not in self._record_parsers:
            raise Exception, "Unrecognized GFF default version: " + version
        self._record_parser = self._record_parsers[version]

        # Stage next record so that we read all metadata at top of file
        self._next_rec = None
        self._stage_rec()

    def get_version(self):
        """Returns the format version used for parsing."""
        return self._version or self._default_version

    def is_version_parsed(self):
        """Returns True if the format version was detected in the stream."""
        return self._version is not None

    def get_metadata(self):
        """Returns the metadata read as a list."""
        return self._metadata

    def get_comments(self):
        """Returns comments read as a list."""
        return self._comments

    def get_sequence_regions(self):
        """Returns the sequence region metadata read as a list."""
        return self._sequence_regions

    def get_fasta_string(self):
        """Returns a string of FASTA formatted sequences found at the end of
        the GFF file (version 3 only)"""
        return self._fasta_string

    def are_references_resolved(self):
        """Returns True if record references have all been resolved."""
        return self._references_resolved

    def read(self):
        """Returns the next record or None if there are none left."""
        try:
            return self.next()
        except StopIteration:
            return None

    def read_recs(self):
        """Returns a list of all records that have not yet been read."""
        return [rec for rec in self]

    def __iter__(self):
        return self

    def next(self):
        self._stage_rec()
        if self._next_rec is None:
            raise StopIteration
        else:
            rec = self._next_rec
            self._next_rec = None
            return rec

    def _stage_rec(self):
        while self._next_rec is None:
            line = self._stream.readline()

            # Stop when EOF reached
            if line == "":
                return
            # Check for pragma line
            elif line.startswith("##"):
                self._parse_directive(line)
            # Check for comment line
            elif line.startswith("#"):
                self._parse_comment(line)
            # Skip over blank lines
            elif line == "\n":
                pass
            # Check for beginning of FASTA region for v3 formats
            elif line.startswith(">") and self._version == "3":
                self._fasta_string = line + self._stream.read()
            else:
                self._next_rec = self._record_parser(line)

    def _set_record_parser(self):
        try:
            self._record_parser = self._record_parsers[self._version]
        except KeyError:
            self._record_parser = self._record_parsers[self._default_version]
            print >>sys.stderr, "Warning: Unrecognized GFF version (%s). " + \
                "Using default version %s." % (self._version, self._default_version)

    def _parse_directive(self, line):
        tokens = line[2:-1].split(None, 1)

        # Skip over empty directive lines
        if not tokens:
            return
            
        # Switch on directive type
        if tokens[0] == "gff-version":
            try:
                self._version = tokens[1]
                self._set_record_parser()
            except IndexError:
                raise FormatError("Invalid gff-version directive: " + line)
        elif tokens[0] == "FASTA":
            # The parser automatically enters FASTA mode with a line starting
            # with '>', so we don't need to do anything here
            pass 
        elif tokens[0] == "#":
            self._references_resolved = True
        # Otherwise this is a metadatum directive
        elif tokens[0] == "sequence-region":
            try:
                seqid, start, end = tokens[1].split()
                seq_region = SequenceRegion(seqid, int(start), int(end))
                self._metadata.append(seq_region)
                self._sequence_regions.append(seq_region)
            except (IndexError, ValueError):
                raise FormatError("Invalid sequence-region directive: " + line)
        else:
            self._metadata.append(Metadatum(*tokens))

    def _parse_comment(self, line):
        # Add line stripped of # prefix and newline to comments list
        self._comments.append(line[1:-1])

    def _parse_record_v1(self, line):
        fields = line[:-1].split('\t', 8)

        if len(fields) == 8:
            attributes = {}
        elif len(fields) == 9:
            attributes = {"group": [fields[8]]}
        else:
            raise FormatError, "Invalid number of fields (should be 8 or 9)"

        try:
            return Record(seqid=fields[0],
                          source=fields[1],
                          type=fields[2],
                          start=int(fields[3]),
                          end=int(fields[4]),
                          score=float(fields[5]),
                          strand=parse_maybe_empty(fields[6]),
                          phase=parse_maybe_empty(fields[7], int),
                          attributes=attributes)
        except ValueError, e:
            raise FormatError, "GFF field format error: " + e.message

    def _parse_record_v2(self, line):
        fields = line[:-1].split('\t', 8)

        if len(fields) == 8:
            attributes_string = ""
        elif len(fields) == 9:
            attributes_string = fields[8]
        else:
            raise FormatError, "Invalid number of fields (should be 8 or 9)"

        try:
            return Record(seqid=fields[0],
                          source=fields[1],
                          type=fields[2],
                          start=int(fields[3]),
                          end=int(fields[4]),
                          score=parse_maybe_empty(fields[5], float),
                          strand=parse_maybe_empty(fields[6]),
                          phase=parse_maybe_empty(fields[7], int),
                          attributes=self._parse_attributes_v2(attributes_string))
        except ValueError, e:
            raise FormatError, "GFF field format error: " + e.message

    def _parse_record_v3(self, line):
        self._references_resolved = False

        fields = line[:-1].split('\t')

        if len(fields) != 9:
            raise FormatError, "Invalid number of fields (should be 9)"

        try:
            return Record(seqid=url_unquote(fields[0]),
                          source=url_unquote(fields[1]),
                          type=url_unquote(fields[2]),
                          start=int(fields[3]),
                          end=int(fields[4]),
                          score=parse_maybe_empty(fields[5], float),
                          strand=parse_maybe_empty(fields[6]),
                          phase=parse_maybe_empty(fields[7], int),
                          attributes=self._parse_attributes_v3(fields[8]))
        except ValueError, e:
            raise FormatError, "GFF field format error: " + e.message

    def _parse_attributes_v3(self, s):
        attributes = {}

        for pair_string in s.split(";"):
            try:
                tag, value = pair_string.split("=")
                attributes[url_unquote(tag)] = map(url_unquote,
                                                   value.split(","))
            except ValueError:
                raise FormatError("Invalid attributes string: " + s)
        return attributes

    def _parse_attributes_v2(self, s):
        return parse_attributes_v2_faster(s)

def parse_attributes_v2(s):
    attributes = {}
    currentTag = None
    for token in AttributeIterator(s):
        if currentTag is None:
            if isinstance(token, IdentifierToken):
                currentTag = token.value
                attributes[currentTag] = []
            else:
                raise FormatError, "Invalid attributes string: " + s
        elif isinstance(token, SeparatorToken):
            currentTag = None
        elif isinstance(token, CommentToken):
            break
        elif isinstance(token, IdentifierToken):
            attributes[currentTag].append(token.value)
        elif isinstance(token, ValueToken):
            attributes[currentTag].append(token.value)
        else:
            raise FormatError, "Invalid attributes string: " + s
    return attributes
    
class IdentifierToken:
    pass
class ValueToken:
    pass
class CommentToken:
    pass
class SeparatorToken:
    pass
class UnknownToken:
    pass

class AttributeIterator:
    identifierPat = re.compile(r'\s*([A-Za-z][A-Za-z0-9_]*)')
    freeTextPat = re.compile(r'\s*"(([^"]|(\\"))*)(?<!\\)"')
    valuePat = re.compile(r'\s*([^;# \t\n\r\f\v]+)')
    sepPat = re.compile(r'\s*(;)')
    commentPat = re.compile(r'\s*#(.*)$')

    pats = (identifierPat, freeTextPat, valuePat, sepPat, commentPat)
    tokenClasses = (IdentifierToken, ValueToken, ValueToken,
                    SeparatorToken, CommentToken)

    def __init__(self, s):
        self.s = s.rstrip()
        self.pos = 0
        
    def __iter__(self):
        return self
    
    def next(self):
        if self.pos >= len(self.s):
            raise StopIteration
        
        for (pat, tclass) in zip(AttributeIterator.pats,
                                 AttributeIterator.tokenClasses):
            match = pat.match(self.s, self.pos)
            if match is not None:
                self.pos = match.end(0)
                t = tclass()
                t.value = match.group(1)
                return t
        else:
            return UnknownToken()

# Faster code for v2 attributes parsing
tag_pat = re.compile(r'\s*([A-Za-z][A-Za-z0-9_]*)')
value_or_sep_pat = re.compile(r'''(?:\s+(?:
                              "((?:[^"]|(?:\\"))*)(?<!\\)"
                              |
                              ([^;#" \t\n\r\f\v]+)
                              ))
                              |
                              \s*(;)''',
                              re.X)
comment_pat = re.compile(r'\s*(?:#(.*))?$')

def parse_attributes_v2_fast(s):
    pos = 0
    d = defaultdict(list)
    while 1:
        # Attempt to match next tag
        tag_match = tag_pat.match(s, pos)
        if tag_match:
            values = d[tag_match.group(1)]
            pos = tag_match.end()
            while 1:
                # Attempt to match values
                value_or_sep_match = value_or_sep_pat.match(s, pos)
                if value_or_sep_match:
                    groups = value_or_sep_match.groups()
                    pos = value_or_sep_match.end()
                    if groups[2]: # separator matched
                        break
                    else:
                        value = groups[0] if groups[0] is not None else groups[1]
                        values.append(value)
                else:
                    raise FormatError, "Invalid attributes string: " + s
        else:
            if not comment_pat.match(s, pos):
                raise FormatError, "Invalid attributes string: " + s
            break
    return d

values_pat = re.compile(r'\s+(?:"((?:[^"]|(?:\\"))*)(?<!\\)")|([^;#" \t\n\r\f\v]+)')
tag_values_pat = re.compile(r'''
                            \s* # optional leading whitespace
                            ([A-Za-z][A-Za-z0-9_]*) # tag
                            (?: # optional value(s)
                            (?: # single value
                            \s+(?:
                            (?:"((?:[^"]|(?:\\"))*)(?<!\\)") # free text value
                            |
                            ([^;#"\s]+) # non-free text value
                            ))
                            |
                            (# multiple values
                            (?:\s+(?:
                            (?:"(?:[^"]|(?:\\"))*(?<!\\)") # free text value
                            |
                            [^;#"\s]+ # non-free text value
                            ))+
                            )
                            )
                            (?:\s*;\s*(?:\#.*$)?) # separator and comment
                            ''',
                            re.X)

def parse_attributes_v2_faster(s):
    tokens = tag_values_pat.split(s)
    d = {}
    for index in xrange(0, len(tokens) - 1, 5):
        if tokens[index]:
            raise FormatError, "Invalid attributes string: " + s
        values = d.setdefault(tokens[index + 1], [])
        if tokens[index + 2] is not None:
            values.append(tokens[index + 2])
        elif tokens[index + 3] is not None:
            values.append(tokens[index + 3])
        else:
            values.extend([q if q else u
                           for q,u in values_pat.findall(tokens[index + 4])])
    if tokens[-1]:
        raise FormatError, "Invalid attributes string: " + s
    return d

def is_integer(x):
    """Returns true if x is of integer type (int or long)."""
    return type(x) in (int, long)

def parse_maybe_empty(s, parse_type=str):
    if s == '.':
        return None
    else:
        return parse_type(s)

def format_maybe_empty(value, empty_str='.'):
    if value is None or value == "":
        return empty_str
    else:
        return str(value)

def quote(s):
    return '"%s"' % str(s)

def url_quote_sub(m):
    return url_quote(m.group(0))

_seqid_pat = re.compile(r'[^a-zA-Z0-9.:^*$@!+_?-|]')
_source_pat = re.compile(r'[^a-zA-Z0-9.: ^*$@!+_?-]')
_type_pat = _source_pat
_tag_pat = re.compile(r'[\t\n\r\f\v;=%&,]')
_value_pat = _tag_pat

class Writer:
    """Writes a GFF formatted file"""

    def __init__(self, stream, version="2", metadata=[]):
        self.stream = stream

        # Create record writer table
        self._record_writers = {"1": self._write_rec_v1,
                                "2": self._write_rec_v2,
                                "2.1": self._write_rec_gtf,
                                "2.2": self._write_rec_gtf,
                                "2.5": self._write_rec_gtf,
                                "3": self._write_rec_v3}

        # Set version
        try:
            self._record_writer = self._record_writers[version]
            self.version = version
            self.write_metadatum(Metadatum("gff-version", version))
        except KeyError:
            raise Exception, "Unrecognized GFF version: " + version

        for metadatum in metadata:
            self.write_metadatum(metadatum)
    
    def write_metadatum(self, metadatum):
        """Writes a metadatum line."""
        if metadatum.value is not None:
            print >>self.stream, "##%s %s" % (metadatum.name, metadatum.value)
        else:
            print >>self.stream, "##%s" % metadatum.name

    def write_comment(self, comment):
        """Writes a comment line."""
        print >>self.stream, "#%s" % comment

    def write(self, rec):
        """Writes a single record."""
        self._record_writer(rec)

    def write_recs(self, recs):
        """Writes a list of records."""
        for rec in recs:
            self.write(rec)

    def _write_rec_v1(self, rec):
        fields = [rec.seqid,
                  rec.source,
                  rec.type,
                  str(rec.start),
                  str(rec.end),
                  format_maybe_empty(rec.score, '0'),
                  format_maybe_empty(rec.strand),
                  format_maybe_empty(rec.phase)]
        if rec.attributes.get('group') is not None:
            fields.append(rec.attributes['group'][0])
        print >>self.stream, '\t'.join(fields)

    def _write_rec_v2(self, rec, attribute_order=None):
        fields = [rec.seqid,
                  rec.source,
                  rec.type,
                  str(rec.start),
                  str(rec.end),
                  format_maybe_empty(rec.score),
                  format_maybe_empty(rec.strand),
                  format_maybe_empty(rec.phase)]
        if rec.attributes:
            if attribute_order is None:
                attribute_order = sorted(rec.attributes.keys())
            fields.append(self._format_attributes_v2(rec.attributes,
                                                     attribute_order))
        print >>self.stream, '\t'.join(fields)

    def _write_rec_v3(self, rec):
        fields = [_seqid_pat.sub(url_quote, rec.seqid),
                  _source_pat.sub(url_quote, rec.source),
                  _type_pat.sub(url_quote, rec.type),
                  str(rec.start),
                  str(rec.end),
                  format_maybe_empty(rec.score),
                  format_maybe_empty(rec.strand),
                  format_maybe_empty(rec.phase),
                  format_maybe_empty(self._format_attributes_v3(rec.attributes))]
        print >>self.stream, '\t'.join(fields)

    def _write_rec_gtf(self, rec):
        # The required GTF attributes
        gtf_attributes = ["gene_id", "transcript_id"]
        
        # Make sure that this record has the required GTF attributes
        if not all([attr in rec.attributes for attr in gtf_attributes]):
            gtf_rec = rec.copy()
            for attr in gtf_attributes:
                gtf_rec.attributes.setdefault(attr, [""])
        else:
            gtf_rec = rec

        attribute_order = gtf_attributes + list(set(gtf_rec.attributes) -
                                                set(gtf_attributes))
        # GTF is just GFF v2 with some required attributes that go in
        # a specific order
        self._write_rec_v2(gtf_rec, attribute_order)
        
    def _format_attributes_v3(self, attributes):
        return ';'.join(["%s=%s" % (_tag_pat.sub(url_quote_sub, tag),
                                    ','.join([_value_pat.sub(url_quote_sub, value)
                                              for value in values]))
                         for tag, values in attributes.items()])

    def _format_attribute_v2(self, tag, values):
        return ' '.join([tag] + map(quote, values)) + ";"

    def _format_attributes_v2(self, attributes, attribute_order):
        return ' '.join([self._format_attribute_v2(tag, attributes[tag])
                         for tag in attribute_order])
