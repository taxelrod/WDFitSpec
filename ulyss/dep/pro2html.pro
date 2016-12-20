;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    PRO2HTML
;       
; PURPOSE:
;    Convert an IDL .pro file into a pretty HTML document, with EMACS-like
;    color coding.
;
; CALLING SEQUENCE:
;    pro2html, pro_file, html_file, signature=
;
; INPUTS: 
;    pro_file: the full path of the IDL .pro file
;    html_file: the full path name of the output html file
;
; OPTIONAL INPUTS:
;    signature: Some text to put at the end.
;
; KEYWORDS:
;    STYLE_SHEET: If present, a string giving the name of an external 
;         style sheet (CSS) that will be linked from the HEAD element of the 
;         generated HTML document. If omitted, a default style will be 
;         included in this document.
;         PRO2HTML defines CSS classes: comment, constant, function-name,
;         keyword, string, variable-name. They can be used to properly
;         colorize the corresponding elements.
;    XREF: If present, a structure used to generate cross-references through
;         the documentation of a package. See ULY_HTML_HELP.
;    PROMPT: If set to a character string, prefix each line with this
;         string as a 'prompt'. 
;    /INCLUDE: If given, only the html body content will be generated
;         (not the page structure) and the extension html_include will be
;         used instead of html.
;         The purpose is to format segment of codes for inclusion in
;         an html document (using SSI for example)
;
; PROCEDURE: 
;    Look for reserved IDL words, comments, and strings and color-code them,
;    converting html senstive characters into codes.  You can change the color
;    coding by changing the style definitions in the header.
;	
;
; REVISION HISTORY:
;    Created: 18-NOV-2000 Erin Scott Sheldon UofMich
;    Inspired by the emacs htmlize.el, moved over to using styles
;    for a more compact and stable result.  E.S. NYU  2006-July-15 
;       
;    2008/07/25, Philippe Prugniel, Observatoire de Lyon:
;    o Replace the call to the procedure dirsep (not found anywhere.
;    o Add the keyword STYLE_SHEET
;-
;
;  This program is part of sdssidl.
;  Copyright (C) 2005  Erin Sheldon, NYU.  erin.sheldon at gmail.com
;
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program; if not, write to the Free Software
;    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
;
;
;                                      



function pro2html_reserved_words

  return,     ['AND', $
               'BEGIN', $
               'BREAK', $
               'CASE', $
               'COMMON', $
               'CONTINUE', $
               'DO',$
               'ELSE', $
               'ELSE:', $
               'END', $
               'ENDCASE',$
               'ENDELSE', $
               'ENDFOR', $
               'ENDIF', $
               'ENDREP', $
               'ENDSWITCH',$
               'ENDWHILE',$
               'EQ', $
               'FOR', $
               'FUNCTION', $
               'GE', $
               'GOTO',$
               'GT', $
               'IF', $
               'INHERITS', $
               'LE', $
               'LT', $
               'MOD', $
               'NE', $
               'NOT', $
               'OF', $
               'ON_ERROR', $
               'ON_IOERROR', $
               'OR', $
               'PRO', $
               'REPEAT', $
               'SWITCH', $
               'THEN',$
               'UNTIL', $
               'WHILE',$
               'XOR','RETURN']


end 

function pro2html_reserved_char

  return, [ ['<', '&lt;'], $
            ['>', '&gt;'], $
            ['&', '&amp;'] ]

end

function is_resword, word

  common pro2html_reser, reswords, html_reschar
  w=where(reswords EQ strupcase(word), nw)
  if nw eq 0 then return, 0 else return, 1

end 


function is_html_reschar, char

  common pro2html_reser, reswords, html_reschar
  w=where(html_reschar[0,*] eq char, nw)
  if nw ne 0 then return, html_reschar[1,w[0]] else return,char

end


pro pro2html, proname, htmlname, silent=silent, signature=signature, $
              style_sheet=style_sheet, xref=xref, prompt=prompt, $
              include=include, header=header

  if n_params() lt 2 then begin 
      print,'-syntax: pro2html, proname, htmlname, silent=silent'
      return
  endif 


  common pro2html_reser, reswords, html_reschar


  if not keyword_set(silent) then begin 
      print
      print,'proname = ',proname
      print,'htmlname = ',htmlname
      print
  endif 


;  dirsep, proname, dir, proc
  dirpro =  strsplit(proname, path_sep(), /EXTR)
  proc = dirpro[n_elements(dirpro)-1]

  openr, inlun, proname, /get_lun

  if keyword_set(include) then exten = '.html_include' $
  else exten = '.html'

  if strtrim(htmlname) eq '' then begin
      file = file_basename(proname,'.pro')+exten
      openw, outlun, file, /get_lun     
  endif else if file_test(htmlname, /DIR) eq 1 then begin
      file = htmlname+path_sep()+file_basename(proname,'.pro')+exten
      openw, outlun, file, /get_lun
  endif else openw, outlun, htmlname, /get_lun

  if n_elements(xref) gt 0 then begin
      if size(xref, /TYPE) eq 8 then begin
          xrf = xref
      endif else begin
          openr, lun, xref, /GET_LUN
          nlines = file_lines(xref)
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
              
              if j eq 0 then xrf = {file:fil, name:sub, purpose:pur}$
              else xrf = [xrf, {file:fil, name:sub, purpose:pur}]
          endfor
          close, lun
          free_lun, lun 
      endelse
  endif

  if not keyword_set(include) then begin
  ;; The header
      printf, outlun, '<html>'
      printf, outlun, '  <head>'
      printf, outlun, '    <title>'+proc+'</title>'
      if n_elements(style_sheet) eq 0 then begin
          printf, outlun, '    <style type="text/css">'
          printf, outlun, '    <!--'
          printf, outlun, '    body {'
          printf, outlun, '      color: #000000;'
          printf, outlun, '      background-color: #ffffff;'
          printf, outlun, '    }'
          printf, outlun, '    .comment {'
          printf, outlun, '      /* font-lock-comment-face */'
          printf, outlun, '      color: #b22222;'
          printf, outlun, '    }'
          printf, outlun, '    .constant {'
          printf, outlun, '      /* font-lock-constant-face */'
          printf, outlun, '      color: #5f9ea0;'
          printf, outlun, '    }'
          printf, outlun, '    .function-name {'
          printf, outlun, '      /* font-lock-function-name-face */'
          printf, outlun, '      color: #0000ff;'
          printf, outlun, '    }'
          printf, outlun, '    .keyword {'
          printf, outlun, '      /* font-lock-keyword-face */'
          printf, outlun, '      color: #a020f0;'
          printf, outlun, '    }'
          printf, outlun, '    .string {'
          printf, outlun, '      /* font-lock-string-face */'
          printf, outlun, '      color: #bc8f8f;'
          printf, outlun, '    }'
          printf, outlun, '    .variable-name {'
          printf, outlun, '      /* font-lock-variable-name-face */'
          printf, outlun, '      color: #b8860b;'
          printf, outlun, '    }'
          printf, outlun
          printf, outlun, '    a {'
          printf, outlun, '      color: inherit;'
          printf, outlun, '      background-color: inherit;'
          printf, outlun, '      font: inherit;'
          printf, outlun, '      text-decoration: inherit;'
          printf, outlun, '    }'
          printf, outlun, '    a:hover {'
          printf, outlun, '      text-decoration: underline;'
          printf, outlun, '    }'
          printf, outlun, '    -->'
          printf, outlun, '    </style>'
      endif else begin
          printf, outlun,'<link rel="stylesheet" type="text/css" href="'+ $
            style_sheet+'">'
       endelse
      Printf,outlun,"<script> pagetitle = '"+strtrim(proc,2)+"' </script>"
      printf, outlun, '  </head><body>'
      printf, outlun, '<script language="JavaScript" type="text/javascript" src="uly.js"></script>'
      if n_elements(header) gt 0 then printf, outlun, header
;      printf,outlun, '<H1 id="page_title">'+proc+'</H1>'
      printf, outlun, '<div style="clear:both;margin:0;padding:0;"></div>'
  endif
  printf, outlun, '<div class="source"><pre>'


  reswords = pro2html_reserved_words()
  html_reschar = pro2html_reserved_char()
  ;; color reserved words
  resf = '<span class="keyword">'
  resb = '</span>'

  ;; color after PRO/FUNCTION statement 
  aprof = '<span class="function-name">'
  aprob = '</span>'

  ;; color comments
  comm = ';'
  comf = '<span class="comment">'
  comb = '</span>'


  ;; color strings dark salmon
  str1 = "'"
  str2 = '"'

  strf = '<span class="string">'
  strb = '</span>'


  nl = file_lines(proname)
  string=''
  for i=0l, nl-1 do begin 

      string=''
      readf, inlun, string
      string=detabify(string)
      length = strlen(string)
      
      if n_elements(prompt) ne 0 then printf, outlun, prompt+' ', F='(A,$)'

      line = ''
      j=0L
      ;; on the line level
      afterpro = 0
      iscomment = 0
      didcomment = 0
      isstring = 0
      oldstr = ''
      while j le length-1 do begin 

          ;; on the word level
          word = ''
          tmp=''
          repeat begin 
              word = word + tmp
              wupcase = strupcase(word)
              if (wupcase eq 'GOTO,') then word=resf+'GOTO'+resb+','
              if (wupcase eq 'RETURN,') then word=resf+'return'+resb+','
              tmp = strmid(string,j,1)
              if afterpro and tmp eq ',' then begin 
                  tmp = aprob+tmp
                  afterpro=0
              endif else if (tmp eq comm) and (not isstring) then begin
                  iscomment = 1 
                  if not didcomment then begin 
                      tmp = comf+tmp
                      didcomment=1
                  endif 
              endif else if ( (tmp eq str1) or (tmp eq str2) ) and (not iscomment) then begin
                  
                  if isstring and (tmp eq oldstr) then begin 
                      isstring = 0 
                      oldstr=''
                      tmp = tmp+strb
                  endif else if not isstring then begin
                      isstring = 1
                      oldstr=tmp
                      tmp = strf+tmp
                  endif 
              endif else tmp = is_html_reschar(tmp)
                  
              j=j+1
          endrep until ( (tmp eq ' ') or (j ge length+1) )

          if iscomment then begin 
              ;; do nothing
          endif else if isstring then begin 
              ;; do nothing
          endif else if afterpro and (tmp ne ' ') then begin 
              word = word + aprob
              afterpro = 0
          endif else if is_resword(word) then begin
              if (strupcase(word) eq 'PRO') or (strupcase(word) eq 'FUNCTION') then afterpro = 1
              if afterpro then word = resf+word+resb+aprof $
              else word = resf+word+resb
          endif 
          line = line + word
          if tmp eq ' ' then line = line + ' '

      endwhile 
      if iscomment then line = line + comb

      if n_elements(xrf) gt 0 then begin ; make the cross-references 
          pos = 0
          while pos ne -1 do begin
              pos = -1
              for kx=0, n_elements(xrf)-1 do begin
                  posi =  stregex(line, xrf[kx].name+'[.,;)( ]', LENGTH=leni)
                  if posi eq -1 then begin  ; search also low-case match
                      posi =  stregex(line, strlowcase(xrf[kx].name)+'[.,;)( ]', LENGTH=leni)
                  endif
                  if posi eq -1 then begin
                      posi =  stregex(line, xrf[kx].name+'$', LENGTH=leni)
                      leni++ ; because this match is shorter than the previous
                  endif
                  if posi eq -1 then begin
                      posi =  stregex(line, strlowcase(xrf[kx].name)+'$', LENGTH=leni)
                      leni++ ; because this match is shorter than the previous
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
                  extr = strmid(line,pos,len-1)
                  printf,outlun, F='(A,$)', strmid(line,0,pos)+$
                    '<a href="'+xrf[km].file+'" onMouseOver="pl('''+xrf[km].purpose+''',+20)" onMouseOut="hl(+20)">'+extr+'</a>'
                  line = strmid(line,pos+len-1)
              endif
          endwhile
      endif

      printf,outlun,line

  endfor 

  printf, outlun, '</pre></div>'
  if not keyword_set(include) then begin
      if n_elements(signature) ne 0 then begin 
          if strtrim(signature) ne '' then printf, outlun, signature
      endif else begin
          printf, outlun, '<hr size=1 noshade>'
          printf, outlun, 'This file was generated from IDL code by pro2html.pro in the SDSSIDL codebase<br>'      
      endelse
      
      printf, outlun, '</BODY>'
      printf, outlun, '</HTML>'
  endif

  free_lun, inlun
  free_lun, outlun

  return
end 
