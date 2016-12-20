;+
; NAME:
;	          ULY_HTML_HELP
;
; PURPOSE:
;                 Generate HTML documentation from the embedded doc
;
; USAGE:
;	          ULY_HTML_HELP, <sources>, <outfile>
;                                [, VERBOSE=<verbose>][, TITLE=<title>]
;                                [, HEADER=<header>]        
;                                [, STYLE_SHEET=<style_sheet>]  
;                                [, CROSSREF=<crossref>]
;                                [, XREF_FILE=<xref_file>]
;
; ARGUMENTS:
;   SOURCES:      A string or string array containing the name(s) of the
;                 .pro files (or the names of directories containing such files)
;                 for which help is desired. If a source file is an IDL procedure,
;                 it must include the .PRO file extension. All other source files
;                 are assumed to be directories.
;
;   OUTFILE:      The name of the output file which will be generated.
;
; KEYWORDS:
;   VERBOSE:	  Normally, ULY_HTML_HELP does its work silently.
;                 Setting this keyword to a non-zero value causes the procedure
;                 to issue informational messages that indicate what it
;                 is currently doing. !QUIET must be 0 for these messages
;                 to appear.
;
;   TITLE:    	  If present, a string which supplies the name that
;                 should appear as the Document and page Title.
;                 By default, the title is "Extended IDL help".
;
;   HEADER:       If present, this string or string-array will be
;                 included at the beginning of the body element,
;                 before the title heading.
;
;   STYLE_SHEET:  If present, a string giving the name of an external 
;                 style sheet (CSS) that will be linked from the HEAD
;                 element of the generated HTML document.
;
;   CROSSREF:     If present, the program will only find the cross-references
;                 and append them to the file specified in XREF_FILE
;
;   XREF_FILE:    If given in combination with CROSSREF, name of the output
;                 file where the cross-reference are appended. Otherwise if
;                 CROSSREF is not set, name of the file used to read the 
;                 cross-references. Cross-references will make the name of 
;                 other documentation blocks clickable.
;                 To use the documentation with cross-references, the
;                 file 'uly.js' must be copied on the same directory as the
;                 html files in order to have the purpose of a function displayed
;                 when the cursor move over its name (this is a
;                 javascript code)
;    
;   KEEPPRO       refs to routine_name.pro not to routine_name.html
;
; DESCRIPTION:
;     This routine is a modified version of the original IDL function
;     MK_HTML_HELP:
;     The differences are:
;        - Make reference to a style-sheet (CSS) that can be passed
;          through the keyword STYLE_SHEET.
;          Define the following classes for some components of the
;          created page:
;                list : list of routines
;                routine : one DIV for each routine documentation 
;          In addition, some elements have particular IDs:
;                page_title : The main title (H1 header)
;                doc_tag : For each tag in the documentation
;        - The list of routines is presented as a table with their
;          name and purpose.
;        - Generate cross references through the different documented modules.
;        - Generate HTMLized files with the source of the program, using
;          the function PRO2HTML.
;
;     Note that PURPOSE is meant to be a short (one line) description
;     of the function. The detailed description should come in the
;     tag DECRIPTION.
;
;     See another variant of the same original program:
;     http://astro.berkeley.edu/~marc/idlshare/general/html/
;
; HISTORY:
;    This routine is a modification of the original IDL Library routine
;    MK_HTML_HELP, Ph. Prugniel, 2008/08

;-
; CATEGORY:       ULY_UTIL
;----------------------------------------------------------------------------
PRO mhh_strict, txtlines
;
; Replaces any occurrence of HTML reserved characters (<,>,&,") in the
; given text lines with the appropriate HTML counterpart.
;
; entry:
;       txtlines - String array containing the text line(s) to be altered.
; exit:
;	txtlines - Same as input except that reserved characters have been 
;                  replaced with the appropriate HTML syntax.
;
 COMPILE_OPT hidden

 count = N_ELEMENTS(txtlines)
 FOR i=0,count-1 DO BEGIN
  txt = txtlines[i] 

  ; Ampersands get replaced with &amp.  Must do ampersands first because
  ; they are used to replace other reserved characters in HTML.
  spos = STRPOS(txt,'&')
  WHILE (spos NE -1) DO BEGIN
   newtxt = STRMID(txt,0,spos)+'&amp;'+STRMID(txt,spos+1,STRLEN(txt)-spos+1)
   txt = newtxt
   spos = STRPOS(txt,'&',spos+1)
  ENDWHILE
  txtlines[i] = txt

  ; Left angle brackets get replaced with &lt;
  spos = STRPOS(txt,'<')
  WHILE (spos NE -1) DO BEGIN
   newtxt = STRMID(txt,0,spos)+'&lt;'+STRMID(txt,spos+1,STRLEN(txt)-spos+1)
   txt = newtxt
   spos = STRPOS(txt,'<',spos+1)
  ENDWHILE
  txtlines[i] = txt

  ; Right angle brackets get replaced with &gt;
  spos = STRPOS(txt,'>')
  WHILE (spos NE -1) DO BEGIN
   newtxt = STRMID(txt,0,spos)+'&gt;'+STRMID(txt,spos+1,STRLEN(txt)-spos+1)
   txt = newtxt
   spos = STRPOS(txt,'>',spos+1)
  ENDWHILE
  txtlines[i] = txt

  ; Double quotes get replaced with &quot;
  spos = STRPOS(txt,'"')
  WHILE (spos NE -1) DO BEGIN
   newtxt = STRMID(txt,0,spos)+'&quot;'+STRMID(txt,spos+1,STRLEN(txt)-spos+1)
   txt = newtxt
   spos = STRPOS(txt,'"',spos+1)
  ENDWHILE
  txtlines[i] = txt

  ; Transform http://... into links 
  spos = stregex(txt,"http://[^ ),]*", LEN=len)
  newtxt = ''
  WHILE (spos NE -1) DO BEGIN
   url = strmid(txt, spos, len)
   newtxt += STRMID(txt,0,spos)+'<a href="'+url+'">'+url+'</a>'
   txt = STRMID(txt,spos+len)
   spos = stregex(txt,"http://[^ ),]*", LEN=len)
  ENDWHILE
  txtlines[i] = newtxt + txt
  
 ENDFOR
END

;----------------------------------------------------------------------------
PRO  mhh_grab_hdr,name,dict,infile_indx,txt_file,verbose,strict
;
; Searches an input file for all text between the ;+ and ;- comments, and
; updates the scratch text file appropriately. Note that this routine
; will extract multiple comment blocks from a single source file if they are
; present.
;
; ARGUMENTS:
;   name          Name of file containing documentation header(s).
;
;   dict[]        Dictionary entries for each documentation block in the .PRO
;                 file.  Each dictionary entry is a structure with an index to 
;                 the source filename, a subject name, scratch file offset,
;                 unique id (for duplicate names), and number of lines of
;                 documentation text.  
;                 This parameter may be originally undefined at entry.
;
;   infile_indx   Index of the source .pro filename.
;
;   txt_file      Scratch file to which the documentation header will
;                 be written.
;
;   verbose       TRUE if the routine should output a descriptive message
;                 when it finds the documentation header.
;
;   strict        If nonzero, the routine will adhere strictly to HTML format.
;                 The document headers will be scanned for characters that are
;                 reserved in HTML (<,>,&,"), which are then converted to the 
;                 appropriate HTML syntax in the output file.
;
; OUTPUTS:
;   txt_file      Updated as necessary. Positioned at EOF.
;
;   dict[]        Updated array of new dictionary entries.
;

 COMPILE_OPT hidden

 ; Under DOS, formatted output ends up with a carriage return linefeed
 ; pair at the end of every record. The resulting file would not be
 ; compatible with Unix. Therefore, we use unformatted output, and
 ; explicity add the linefeed, which has ASCII value 10.
 LF=10B

 OPENR, in_file, /GET, name

 IF (verbose NE 0) THEN MESSAGE,/INFO, 'File = '+name
 WHILE (1) DO BEGIN
  ; Find the opening line of the next header.
  tmp = ''
  found = 0
  num = 0
  header = ''
  ON_IOERROR, DONE
  WHILE (NOT found) DO BEGIN
   READF, in_file, tmp
   IF (STRMID(tmp,0,2) EQ ';+') THEN found = 1
  ENDWHILE

  IF (found) THEN BEGIN
   ; Find the matching closing line of the header.
   found = 0
   WHILE (NOT found) DO BEGIN
    READF,in_file,tmp
    IF (STRMID(tmp,0,2) EQ ';-') THEN BEGIN
     found =1
    ENDIF ELSE BEGIN
     tmp = strmid(tmp, 1, 1000)
     header = [header, tmp]
     num = num + 1
    ENDELSE
   ENDWHILE

   IF (strict) THEN mhh_strict,header
   ; Done with one block of header

   ; Keep track of current scratch file offset, then write doc text.
   POINT_LUN,-txt_file,pos
   FOR i=1, num DO BEGIN
    WRITEU, txt_file, header[i],LF
   ENDFOR

   ; Search for the subject. It is the line following name.
   ; Alteration for ULySS below ... take the 'subject' even if it is
   ; on the same line as 'NAME:'
   index = WHERE(STRTRIM(header, 2) EQ 'NAME:', count)
   IF (count eq 1) THEN BEGIN
    sub = STRUPCASE(STRTRIM(header[index[0]+1], 2))
    IF (verbose NE 0) THEN MESSAGE,/INFO, 'Routine = '+sub
   ENDIF ELSE BEGIN
       index = where(strmatch(strtrim(header,2),'NAME: *'), count)
       if (count eq 1) then begin
           sub = strupcase(strtrim(strmid(header[index[0]],6),2))
           IF (verbose NE 0) THEN MESSAGE,/INFO, 'Routine = '+sub
       endif
   endelse

   ; If the NAME field was not present, set the subject to the name of the 
   ; source text file.
   if count ne 1 then begin
    MESSAGE,/INFO,'Properly formatted NAME entry not found...'
    ifname = name

    tok = PATH_SEP()

    ; Cut the path.
    sp0 = 0
    spos = STRPOS(ifname,tok,sp0)
    WHILE (spos NE -1) DO BEGIN
     sp0 = spos+1
     spos = STRPOS(ifname,tok,sp0)
    ENDWHILE
    ifname = STRMID(ifname,sp0,(STRLEN(ifname)-sp0))

    ; Cut the suffix.
    spos = STRPOS(ifname,'.')
    IF (spos NE -1) THEN ifname = STRMID(ifname,0,spos[0])
    IF (strict) THEN mhh_strict, ifname
    sub = STRUPCASE(ifname)
    MESSAGE,/INFO,'  Setting subject to filename: '+sub+'.'
   ENDif

   ; Extract the purpose line (yet assume it is on a single line)
   purpose=''
   index = WHERE(STRTRIM(header, 2) EQ 'PURPOSE:', count)
   IF (count eq 1) THEN BEGIN
    purpose = STRTRIM(header[index[0]+1], 2)
   ENDIF ELSE BEGIN
       index = where(strmatch(strtrim(header,2),'PURPOSE: *'), count)
       if (count eq 1) then begin
           purpose = strtrim(strmid(header[index[0]],9))
       endif
   endelse

   ; Calculate unique id in case of duplicate subject names.
   IF (N_ELEMENTS(dict) EQ 0) THEN $
    ndup=0 $
   ELSE BEGIN
    dpos = WHERE(dict.subject EQ sub,ndup)
    IF (ndup EQ 1) THEN dict[dpos[0]].id = 1
    IF (ndup GE 1) THEN ndup = ndup + 1
   ENDELSE

   ; Create a dictionary entry for the document header.
   entry = {DICT_STR,subject:sub,purpose:purpose,indx:infile_indx, $
            id:ndup,offset:pos,nline:num}
   IF (N_ELEMENTS(dict) EQ 0) THEN dict = [entry] ELSE dict = [dict,entry]
  ENDIF
 ENDWHILE

DONE: 
 ON_IOERROR, NULL
 FREE_LUN, in_file
END

;----------------------------------------------------------------------------
PRO mhh_gen_file,dict,txt_file,infiles,outfile,verbose,title,strict,style_sheet, xref,header, KEEPPRO=keeppro
; DESCRIPTION:
;     Build a .HTML file with the constituent parts.
;
; ARGUMENTS:
;   dict         Array of dictionary entries. Each entry is a structure
;                with a subject name, scratch file offset, number of lines
;                of documentation text, etc.
;
;   infiles      String array containing the name(s) of .pro files 
;                for which help is being generated.
;
;   txt_file     Scratch file containing the documentation text.
;
;   outfile      NAME of final HELP file to be generated.
;
;   verbose      TRUE if the routine should output a descriptive message
;                when it finds the documentation header.
;
;   title        Scalar string containing the name to go at the top of the
;                HTML help page.
;
;   strict       If nonzero, the routine will adhere strictly to HTML format.
;                The document headers will be scanned for characters that are
;                reserved in HTML (<,>,&,"), which are then converted to the 
;                appropriate HTML syntax in the output file.
;
; OUTPUTS:
;   outfile      It has been created.
;
;   txt_file     It has been closed via FREE_LUN.
;

 COMPILE_OPT hidden

 IF (N_ELEMENTS(dict) EQ 0) THEN BEGIN
    MESSAGE, 'No documentation headers found', LEVEL=-1, /CONTINUE
    PRINT, ' '
    PRINT, 'None of the specified files contain documentation header'
    PRINT, 'information in the DOC_LIBRARY format. See the documentation'
    PRINT, 'for ULY_HTML_HELP and DOC_LIBRARY for details.'
    PRINT, ' '
    PRINT, 'No HTML files were created.'
    PRINT, ' '
    IF((N_ELEMENTS(txt_file) NE 0) && ((FSTAT(txt_file)).open EQ 1)) $
       THEN FREE_LUN, txt_file
    RETURN
 ENDIF

 ; Append unique numbers to any duplicate subject names.
 dpos = WHERE(dict.id GT 0,ndup) 
 FOR i=0,ndup-1 DO BEGIN
  entry = dict[dpos[i]]
  dict[dpos[i]].subject = entry.subject+'['+STRTRIM(STRING(entry.id),2)+']'
 ENDFOR

 ; Sort the subjects alphabetically.
 count = N_ELEMENTS(dict)
 indices = SORT(dict.subject)

 ; Open the final file.
 OPENW,final_file,outfile,/GET_LUN
 IF (verbose NE 0) THEN MESSAGE,/INFO,'Building '+outfile+'...'

 ; Print a comment indicating how the file was generated.
 PRINTF,final_file,'<!-- This file was generated by uly_html_help.pro -->'

 ; Header stuff.
 PRINTF,final_file,'<html>'
 PRINTF,final_file,' '

 PRINTF,final_file,'<head>'
 if n_elements(style_sheet) ne 0 then begin                      ; Style sheet
   PRINTF,final_file,'<link rel="stylesheet" type="text/css" href="'+ $
     style_sheet+'">'
endif
 PRINTF,final_file,'<TITLE>',title,'</TITLE>'                        ; Title.
 PRINTF,final_file,"<script> pagetitle = '"+strtrim(title,2)+"' </script>"
 PRINTF,final_file,'</head>'
 PRINTF,final_file,' '

 ; Title and intro info.
 PRINTF,final_file,'<body>'
 if n_elements(header) gt 0 then PRINTF,final_file, header
 PRINTF,final_file,'<script language="JavaScript" type="text/javascript" src="uly.js"></script>'
; PRINTF,final_file,'<H1 id="page_title">',title,'</H1>'
 PRINTF,final_file,'<P style="clear:both" ></P>'

 ; Index.
 PRINTF,final_file,'<A NAME="ROUTINELIST">'
 PRINTF,final_file,'<H1>List of Routines</H1></A>'
 PRINTF,final_file,'<DIV CLASS="list">'
; PRINTF,final_file,'<UL>'
 PRINTF,final_file,'<TABLE>'

 FOR i=0,count-1 DO BEGIN
  entry = dict[indices[i]]

  IF (entry.nline GT 0) THEN $
   PRINTF,final_file,'<tr><td><A HREF="#'+entry.subject+'">',entry.subject,'</A>'+'</td><td>'+entry.purpose+'</td></tr>'

 ENDFOR
; PRINTF,final_file,'</UL><P>'
 PRINTF,final_file,'</TABLE>'
 PRINTF,final_file,' '

 PRINTF,final_file,'</DIV>'
 PRINTF,final_file,'<HR>'
 PRINTF,final_file,' '

 ; Descriptions.
 PRINTF,final_file,'<H1>Routine Descriptions</H1>'
 ON_IOERROR,TXT_DONE
 FOR i=0,count-1 DO BEGIN
  entry = dict[indices[i]]
  IF (entry.nline GT 0) THEN BEGIN
   PRINTF,final_file,'<DIV CLASS="routine">'
   PRINTF,final_file,'<A NAME="',strtrim(entry.subject,2),'">'
   PRINTF,final_file,'<H2>',entry.subject,'</H2></A>'

   prev_i = i - 1
   IF (prev_i LT 0) THEN $
    dostep = 0 $ 
   ELSE BEGIN
    prev_ent = dict[indices[prev_i]]
    dostep = prev_ent.nline EQ 0
   ENDELSE
   WHILE dostep DO BEGIN
    prev_i = prev_i - 1
    IF (prev_i LT 0) THEN $
     dostep = 0 $
    ELSE BEGIN
     prev_ent = dict[indices[prev_i]]
     dostep = prev_ent.nline EQ 0
    ENDELSE
   ENDWHILE
   IF (prev_i GE 0) THEN $
    PRINTF,final_file,'<A HREF="#',prev_ent.subject,'">[Previous Routine]</A>'

   next_i = i + 1
   IF (next_i GE count) THEN $
    dostep = 0 $
   ELSE BEGIN
    next_ent = dict[indices[next_i]]
    dostep = next_ent.nline EQ 0
   ENDELSE
   WHILE dostep DO BEGIN
    next_i = next_i + 1
    IF (next_i GE count) THEN $
     dostep = 0 $
    ELSE BEGIN
     next_ent = dict[indices[next_i]]
     dostep = next_ent.nline EQ 0
    ENDELSE
   ENDWHILE
   IF (next_i LT count) THEN $
    PRINTF,final_file,'<A HREF="#',next_ent.subject,'">[Next Routine]</A>'

   PRINTF,final_file,'<A HREF="#ROUTINELIST">[List of Routines]</A>'
   PRINTF,final_file,'<PRE>'
   tmp = ''

   POINT_LUN,txt_file,entry.offset
   FOR j=1,entry.nline DO BEGIN
    READF,txt_file,tmp
    s = STREGEX(tmp, '^ [1-Z][A-Z0-9_ ]*:',LENGTH=len)   ; search doc tags
    if len gt 0 then begin
      printf,final_file, '<span id="doc_tag">',strmid(tmp,0,len),'</span>',strmid(tmp,len)
    endif else begin                               ; check the xref...
        pos = 0
        while pos ne -1 do begin
            pos = -1
            for kx=0, n_elements(xref)-1 do begin
                posi =  stregex(tmp, xref[kx].name+'[.,;) ]', LENGTH=leni)
                if posi eq -1 then begin
                    posi =  stregex(tmp, xref[kx].name+'$', LENGTH=leni)
                    leni++   ; because this match is shorter than the previous
                endif
                if posi ne -1 then begin
                    if pos eq -1 then pos=posi $
                    else if posi lt pos then pos = posi
                    if pos eq posi then begin
                        km = kx
                        len = leni
                    endif
                endif
            endfor
            if pos ne -1 then begin
                extr = strmid(tmp,pos,len-1)
                PRINTF,final_file, F='(A,$)', strmid(tmp,0,pos)+$
                  '<a href="'+xref[km].file+'" onMouseOver="pl('''+xref[km].purpose+''',+20)" onMouseOut="hl(+20)">'+extr+'</a>'
                tmp = strmid(tmp,pos+len-1)
            endif
        endwhile
        PRINTF,final_file,tmp
    endelse
   ENDFOR
   PRINTF,final_file,'</PRE><P>'

   dirpro = strsplit(outfile, path_sep(), /EXTR)
   if n_elements(dirpro) gt 1 then begin
       dir = strjoin(dirpro[0:n_elements(dirpro)-2], path_sep()) + path_sep()
   endif else dir=''
   dirpro = strsplit(infiles[entry.indx], path_sep(), /EXTR)
   if not keyword_set(keeppro) then begin
      progn = dirpro[n_elements(dirpro)-1]
      progn = strmid(progn,0,strlen(progn)-4) + '.html'
      pro2html, infiles[entry.indx], dir+progn, $
         /SILENT, $
         STYLE_SHEET=style_sheet, $
         XREF=xref, $
         HEADER = header, $
         SIGN='Part of <a href=".">ULySS package</a>'
      fname = infiles[entry.indx]
   endif else begin
      fname = dirpro[n_elements(dirpro)-1]
      progn = fname
   endelse

   IF (strict) THEN mhh_strict,fname

   fname = '<A HREF="' + progn + '">' + fname + '</A>'

   PRINTF,final_file,'<STRONG>(See '+fname+')</STRONG><P>'
   PRINTF,final_file,'</DIV>'
   PRINTF,final_file,'<HR>'
   PRINTF,final_file,' '
  ENDIF
 ENDFOR
TXT_DONE:
 ON_IOERROR,NULL
 FREE_LUN,txt_file

 PRINTF,final_file,'This page was created by the ULySS routine '
 PRINTF,final_file,'<CODE>uly_html_help</CODE>.<BR>'
 PRINTF,final_file,'<TABLE WIDTH="97%"><TR>'
 PRINTF,final_file,'<TD ALIGN="left">'
 PRINTF,final_file,'<STRONG>Contact:</STRONG> <I>ulyss at obs.univ-lyon1.fr</I>'
 PRINTF,final_file,'</TD><TD ALIGN="right">'
 PRINTF,final_file,'<STRONG>Last modified: </STRONG><I>',SYSTIME(),'</I>.<BR>'
 PRINTF,final_file,'</TD></TR></TABLE>'
 PRINTF,final_file,' '
 PRINTF,final_file,'<HR>'
 PRINTF,final_file,' '

 ; Footer.
 PRINTF,final_file,'</body>'
 PRINTF,final_file,'</html>'
 FREE_LUN,final_file

 ; Tell the user the file was created.
 MESSAGE, 'Created file '+(FILE_INFO(outfile)).name, /INFORMATIONAL, LEVEL=-1
END

;----------------------------------------------------------------------------
PRO uly_html_help, sources, outfile, VERBOSE=verbose, TITLE=title, HEADER=header, STYLE_SHEET=style_sheet, CROSSREF=crossref, XREF_FILE=xref_file, KEEPPRO=keeppro

 ON_ERROR, 2

 strict = 1   ; Suppress the non-strict mode of MK_HTML_HELP

 ; Establish a catch block so that we can close the scratch file on error.
 catch, err
 if (err ne 0) then begin
    catch,/cancel
    ; If we got as far as opening the scratch file, then close it.
    ; It will be deleted on close because it was opened with /DELETE.
    IF((n_elements(txt_file) ne 0) && ((FSTAT(txt_file)).open EQ 1)) $
	THEN FREE_LUN,txt_file
    MESSAGE,/REISSUE_LAST
 endif

 IF (NOT KEYWORD_SET(verbose)) THEN verbose=0
 IF (NOT KEYWORD_SET(title)) THEN title="Extended IDL Help" 

 infiles = ''

 count = N_ELEMENTS(sources)
 IF (count EQ 0) THEN BEGIN
  MESSAGE,/INFO,'No source IDL directories found.'
  RETURN
 ENDIF

 ; Open a temporary file for the documentation text.
 OPENW, txt_file, FILEPATH('userhtml.txt', /TMP), /GET_LUN, /DELETE

 ; Loop on sources. 
 FOR i=0, count-1 DO BEGIN
  src = sources[i]

  ; Test if the file is a .PRO file.
  IF (STRUPCASE(STRMID(src, STRLEN(src)-4,4)) EQ '.PRO') THEN BEGIN 
    infiles = [infiles,src]
  ENDIF ELSE BEGIN	; Not a .PRO file, treat it as a directory.
    ; Get a list of all .pro files in the directory.
    flist = FILE_SEARCH(src+PATH_SEP()+'*.pro',COUNT=npro)
    IF (npro GT 0) THEN infiles = [infiles,flist]
  ENDELSE
 ENDFOR

 count = N_ELEMENTS(infiles)
 IF (count EQ 1) THEN BEGIN
  MESSAGE,/INFO,'No IDL files found.'
  IF((FSTAT(txt_file)).open EQ 1) THEN FREE_LUN,txt_file
  RETURN
 ENDIF 
 infiles = infiles[1:*]
 count = count-1

 ; Loop on all files.
 FOR i=0,count-1 DO BEGIN
   name = infiles[i]
   mhh_grab_hdr,name,dict,i,txt_file,verbose,strict
 ENDFOR

; generate cross-references
 if keyword_set(crossref) then begin
     if n_elements(outfile) gt 0 then begin
         outf = strsplit(outfile,PATH_SEP(), /EXTRACT)
         outf = outf[n_elements(outf)-1]
     endif else outf = 'index.html'
     if n_elements(xref_file) gt 0 then $
       openw, lun, xref_file, BUFSIZE=0, /GET_LUN, /APPEND $
     else lun = -1
     for k=0,n_elements(dict)-1 do begin
         printf, lun, outf+'#'+strtrim(dict[k].subject,2)+' '+strtrim(dict[k].subject,2)+' '+dict[k].purpose
     endfor
     if n_elements(xref_file) gt 0 then begin
         close, lun
         free_lun, lun 
     endif
     return
 endif
 
; load the cross-reference file
 if n_elements(xref_file) eq 1 then begin
     openr, lun, xref_file, /GET_LUN
     nlines = file_lines(xref_file)
     strbuf=strarr(1)
     for j = 0L, nlines-1 do begin
         readf, lun, strbuf
         cols = strsplit(strbuf)
         fil = strtrim(strmid(strbuf,0,cols[1]-1),2)
         if n_elements(cols) ge 3 then begin
             sub = strtrim(strmid(strbuf,cols[1],cols[2]-cols[1]-1),2)
             pur = strtrim(strmid(strbuf,cols[2]),2)
         endif else if cols eq 2 then begin
             sub = strtrim(strmid(strbuf,cols[1]),2)
             pur = ''
         endif else message, 'Wrong xref file'
           
         if j eq 0 then xref = {file:fil, name:sub, purpose:pur}$
         else xref = [xref, {file:fil, name:sub, purpose:pur}]
     endfor
     close, lun
     free_lun, lun 
 endif

 ; Generate the HTML file.
 mhh_gen_file,dict,txt_file,infiles,outfile,verbose,title,strict,style_sheet, xref, header, KEEPPRO=keeppro

; dirpro = strsplit(outfile, path_sep(), /EXTR)
; if n_elements(dirpro) gt 1 then begin
;     dir = strjoin(dirpro[0:n_elements(dirpro)-2], path_sep()) + path_sep()
; endif else dir=''
; for k=0,n_elements(infiles)-1 do begin
;     dirpro = strsplit(infiles[k], path_sep(), /EXTR)
;     progn = dirpro[n_elements(dirpro)-1]
;     progn = dir + strmid(progn,0,strlen(progn)-4) + '.html'
;     print, progn
; endfor

END
