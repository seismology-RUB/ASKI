#!/usr/bin/env python
#
#----------------------------------------------------------------------------
#   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
#
#   This file is part of ASKI version 1.2.
#
#   ASKI version 1.2 is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.
#
#   ASKI version 1.2 is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.
#----------------------------------------------------------------------------
#
import argparse
from sys import exit as sys_exit
#
characters_to_check_header_start = '----------------------------------------------------------------------------'
characters_to_check_header_end =   '----------------------------------------------------------------------------'
GPL_header_line_segments = ['under the terms of the GNU General Public License',
                            'You should have received a copy of the GNU General Public License',
                            'under the terms of the GNU Free Documentation License']
#
#--------------------------------------------------------------------------------
# AUXILIARY FUNCTIONS
#--------------------------------------------------------------------------------
def exit_now():
    print "\n+++ Good Bye\n"
    sys_exit()
#
#--------------------------------------------------------------------------------
#
def replace_first_word_after_key(original_string,key,replace_string):
    string_split = original_string.split(key)
    word_after_key = string_split[1].split()[0]                  # first word behind first occurrence of key
    if word_after_key.endswith('.'):
        # if word_after_key concludes a sentence (by period / full stop), do NOT replace the period (only the last period in word_after_key is preserved)
        word_after_key = word_after_key[:-1]
    return (string_split[0]+                                     # everything in front of first occurence of key in original_string
            key+                                                 # this was removed by .split(key)
            string_split[1].split(word_after_key)[0]+            # everything (e.g. whitespace) in front of first occurence of key
            replace_string+                                      # put replace string instead of first occurence of key
            word_after_key.join(string_split[1].split(word_after_key)[1:])+  # exact portion of original string following word_after_key
            key.join(string_split[2:]))                          # in case there are more occurences of key in original_String (i.e. len(string_split)>2), put the rest of original_string
#
#--------------------------------------------------------------------------------
#
def change_header(filename,nlines,year=None,version=None):
    adapt_year = year is not None
    adapt_version = version is not None
    print "file '"+filename+"'"

    # try to read in the complete file content (is required for writing edited file below)
    # IN THE FUTURE: if it is not feasible to read in complete files (e.g. too much memory
    # required), you need to re-implement this routine, only reading in the first nlines
    # and modifying them without reading the rest of the lines in advance (e.g. looping on all
    # remaining lines, etc.)
    try:
        lines = open(filename,'r').readlines()
    except:
        print "   WARNING: could not open file; file is ignored"
        return

    if len(lines) < nlines:
        print "   WARNING: could not read "+str(nlines)+" lines from file; file is ignored"
        return

    # assume GPL header to start with the first line containing 76 '-' characters
    # iline_header_start points to the NEXT line
    iline_header_start = min([i+1 if characters_to_check_header_start in line
                              else nlines for i,line in enumerate(lines[:nlines])])
    # assume GPL header to end with the first line containing 76 '-' characters AFTER the header-starting-line
    # iline_header_end points to the line BEFORE that
    if iline_header_start < nlines:
        iline_header_end = min([i+iline_header_start if characters_to_check_header_end in line
                                else nlines for i,line in enumerate(lines[iline_header_start:nlines])])
        # if there was no end line found (all entries equal nlines), set iline_header_end to a value 
        # that creates an empty loop below
        if iline_header_end == nlines:
            iline_header_end = iline_header_start
    else:
        iline_header_end = iline_header_start

    if iline_header_end <= iline_header_start:
        print "   WARNING: could not find start AND end line within first "+str(nlines)+" lines of the file"
        return

    GPL_header_line_segment_was_found = any([GPL_segment in line for GPL_segment in GPL_header_line_segments
                                             for line in lines[iline_header_start:iline_header_end]])
    if not GPL_header_line_segment_was_found:
        print ("   WARNING: could not find a line segment indicating GPL (or GFDL) between "+
               "start and end line (lines "+str(iline_header_start)+" to "+str(iline_header_end+1)+")file is ignored")
        return

    # after all checks being succesful (GPL header found), now loop on header lines and adapt those as requested
    year_was_changed = False
    version_was_changed = False
    for i,line in enumerate(lines[iline_header_start:iline_header_end]):
        iline = i+iline_header_start
        #print iline,line.strip()

        # the year is assumed to be contained only in first line after start line of header after word 'Copyright '
        if i==0:
            if adapt_year:
                if 'Copyright ' in lines[iline]:
                    try:
                        line_adapted = replace_first_word_after_key(line.strip('\n'),'Copyright ',year)+"\n"
                        #print line_adapted.strip()
                        lines[iline] = line_adapted
                        year_was_changed = True
                    except:
                        print "   WARNING: could not replace next word after first 'Copyright ' in line "+str(iline+1)
        # search the rest of the lines to search for 'ASKI version '
        else:
            if adapt_version:
                if 'ASKI version ' in lines[iline]:
                    try:
                        line_adapted = replace_first_word_after_key(line.strip('\n'),'ASKI version ',version)+"\n"
                        #print line_adapted.strip()
                        lines[iline] = line_adapted
                        version_was_changed = True
                    except:
                        print "   WARNING: could not replace next word after first 'ASKI version ' in line "+str(iline+1)

    if adapt_year and not year_was_changed:
        print "   WARNING: could not adapt year (though requested): keyword 'Copyright' not present in first line of header after start line"
    if adapt_version and not version_was_changed:
        print "   WARNING: could not adapt version (though requested): keyword 'ASKI version' nowhere present in GPL header"

    file_was_changed = year_was_changed or version_was_changed
    if file_was_changed:
        try:
            open(filename,'w').writelines(lines)
        except:
            print "   WARNING: could not overwrite file with adapted file content"

        message = "   the GPL header was successfully adapted (by "
        if year_was_changed:
            message += "year"
        if year_was_changed and version_was_changed:
            message += " and "
        if version_was_changed:
            message += "version"
        message += ")"
        print message
    # else: # file_was_changed
    #     # in this case year_was_changed == False and version_was_changed == False
    #     # so specify warning message solely on what was requested (here this is not fulfilled)
    #     message = "   WARNING: could not adapt the GPL header of the file : "
    #     if adapt_year:
    #             message+="keyword 'Copyright' not present in first line of header after start line; "
    #     if adapt_version:
    #             message+="keyword 'ASKI version' nowhere present in GPL header ; "
    #     print message
#
#--------------------------------------------------------------------------------
# MAIN PROGRAM
#--------------------------------------------------------------------------------
def main():
    # initialize argument parser 
    parser = argparse.ArgumentParser(prog='adapt_GPL_headers.py',
                                     description='for ASKI source files given as positional arguments (at least 1 '+
                                     'must be given), adapt GPL header as specified (year, ASKI version)')
    # add options
    parser.add_argument('-y','--year', required=False, help='if set, then copyright year (second word of GPL header) '+
                        'will be set to given value')
    parser.add_argument('-ver','--version', required=False, help='if set, then ASKI version will be set to given value: '+
                        'ONLY give numbers here like "1.3" for ASKI version 1.3')
    parser.add_argument('-n','--nlines', required=False, default=21, type=int, help='if set, the GPL header is expected to '+
                        'be contained in the first n lines of each file (default: 21)')
    # add positional arguments
    parser.add_argument('ASKI_src_files', nargs='+', help='list of ASKI source files (at least one file must be given)')
    
    args = parser.parse_args()
    
    if not (args.year or args.version):
        parser.error('No action requested, please indicate at least one of options -y or -ver')
        
    if args.year:
        #print args.year
        print "year will be changed in: '"+args.year+"'"
    else:
        print "no year given: year will not be changed"

    if args.version:
        #print args.version
        print "version will be changed in: '"+args.version+"'"
    else:
        print "no version given: version will not be changed"

    print "files in which GPL header will be adapted (if there is a header):\n"+' '.join(args.ASKI_src_files)+"\n"

    print "assuming GPL headers within the first "+str(args.nlines)+" lines of a file\n"

    print "GPL headers are defined by:"
    print ("   start line: first line containing the following characters (form of start line):\n"+
           "      '"+characters_to_check_header_start+"'")
    print ("   end line: first line after starting line containing the following characters (form ofend line):\n"+
           "      '"+characters_to_check_header_start+"'")
    print ("   there is ONE line segment out of the following, somewhere between start and end line:\n"+
           "      '"+"'\n      '".join(GPL_header_line_segments)+"'")
    print ("   Assumption: first line after start line does NOT contain 'ASKI version ' (ignored)\n"+
           "      but contains 'Copyright ' followed by the year (following word will be\n"+
           "      replaced by year), nowhere else will the year be adapted")
    print ""

    for filename in args.ASKI_src_files:
        change_header(filename,args.nlines,args.year,args.version)

    exit_now()


if __name__=="__main__":
    main()
