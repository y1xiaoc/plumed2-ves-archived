" Vim syntax file
" Language: PLUMED

if exists("b:current_syntax")
  finish
endif

let b:current_syntax="plumed"

if exists("g:plumed_shortcuts")
  noremap  <buffer> <F2> :PHelp<CR>
  inoremap <buffer> <F2> <Esc>:PHelp<CR>
endif

" path for plumed plugin
let s:path=expand('<sfile>:p:h:h')

" All except space and hash are in word
setlocal iskeyword=33,34,36-126

" Matching dots, possibly followed by a comment
" Only dots are part of the match
syntax match plumedDots /\v\.\.\.(\s*(#.*)*$)@=/ contained
highlight link plumedDots Type

let b:plumedActions=[]
let b:plumedDictionary={}

function! PlumedDefineSyntax()

  for key in sort(keys(b:plumedDictionary))
    call add(b:plumedActions,{"word":key})
  endfor

for a in b:plumedActions
  let action=a["word"]
" vim variables cannot contain -
" we convert it to triple ___
  let action_=substitute(action,"-","___","g")

  for b in b:plumedDictionary[action]
    if(b["menu"]=="(numbered)")
      let string='"\v<' . b["word"] . '[0-9]*\=[^{ #]*"'
    elseif(b["menu"]=="(option)")
" this is necessary since word for option is e.g ."RESTART="
      let string='"\v<' . substitute(b["word"],"=","","") . '[0-9]*\=[^{ #]*"'
    elseif(b["menu"]=="(flag)")
      let string='"\v<' . b["word"] . '>"'
    endif
    execute 'syntax match   plumedKeywords' . action_ . ' ' . string . ' contained contains=plumedStringInKeyword'
  endfor

" single line, with explicit LABEL
" matching action at beginning of line, till the end of the line
" can contain all the keywords associated with this action, plus strings, label, and comments
execute 'syntax region plumedLine' . action_ . ' matchgroup=plumedAction' . action_ . ' start=/\v^\s*' . action . '>/ excludenl end=/$/ contains=plumedComment,plumedKeywords' . action_ . ',plumedLabel,plumedStringOneline fold'
" multiple line, with explicit LABEL
" first row might contain extra words before arriving at the dots
" thus continuation dots are matched by plumedDots
" matching action at beginning of line, followed by dots and possibly by a comment
" ends on dots, possibly followed by the same action name and possibly a comment
" comments and initial dots are not part of the match
" can contain all the keywords associated with this action, plus strings, label, and comments
execute 'syntax region plumedCLine' . action_ . ' matchgroup=plumedAction' . action_ . ' start=/\v^\s*' . action . '>(.+\.\.\.\s*(#.*)*$)@=/ end=/\v^\s*\.\.\.(\s+' . action . ')?\s*((#.*)*$)@=/ contains=plumedComment,plumedKeywords' . action_ . ',plumedLabel,plumedString,plumedDots fold'
" single line, with label: syntax
" matching label followed by action
" can contain all the keywords associated with this action, plus strings and comments
" labels are not allwed
execute 'syntax region plumedLLine' . action_ . ' matchgroup=plumedAction' . action_ . ' start=/\v^\s*[^ #@][^ #]*:\s+' . action . '/ excludenl end=/$/ contains=plumedComment,plumedKeywords' . action_ . ',plumedStringOneline fold'
" multiple line, with label: syntax
" first row might contain extra words before arriving at the dots
" thus continuation dots are matched by plumedDots
" matching label, action, dots, and possibly comment
" comments and dots are not part of the match
" ends on dots, possibly followed by the same label and possibly a comment
" comments and initial dots are not part of the match
execute 'syntax region plumedLCLine' . action_ . ' matchgroup=plumedAction' . action_ . ' start=/\v^\s*\z([^ #@][^ #]*\:)\s+' . action . '>(.+\.\.\.\s*(#.*)*$)@=/ end=/\v^\s*\.\.\.(\s+\z1)?\s*((#.*)*$)@=/ contains=plumedComment,plumedKeywords' . action_ . ',plumedString,plumedDots fold'
" this is a hack required to match the ACTION when it is in the second line
execute 'syntax match plumedSpecial' . action_ . ' /\v(\.\.\.\s*(#.*)*\_s*)@<=' . action . '>/ contained'
execute 'highlight link plumedSpecial' . action_ . ' Type'
" multiple line, with label: syntax
" here ACTION is on the second line
" matching label, dots, possibly comments, newline, then action name
" comments, dots, and action are not part of the match
" ends on dots possibly followed by the same label and possibly a comment
execute 'syntax region plumedLCLine' . action_ . ' matchgroup=plumedAction' . action_ . ' start=/\v^\s*\z([^ #@][^ #]*\:)\s+(\.\.\.\s*(#.*)*\_s*' . action . ')@=/ end=/\v^\s*\.\.\.(\s+\z1)?\s*((#.*)*$)@=/ contains=plumedComment,plumedKeywords' . action_ . ',plumedString,plumedSpecial' . action_ . ',plumedDots fold'
execute 'highlight link plumedAction' . action_ . ' Type'
execute 'highlight link plumedKeywords' . action_ . ' Statement'
endfor

" comments and strings last, with highest priority
syntax region  plumedString start=/\v\{/  end=/\v\}/ contained contains=ALL fold
syntax region  plumedStringOneline start=/\v\{/  end=/\v\}/ oneline contained contains=ALL fold
highlight link plumedString String
highlight link plumedStringOneline String
syntax match   plumedStringInKeyword /\v(<[^ #]+\=)@<=[^ #]+/ contained
highlight link plumedStringInKeyword String

" Matching label
syntax match   plumedLabel "\v<LABEL\=[^ #]*" contained contains=plumedLabelWrong
highlight link plumedLabel Type

" Errors
syntax match   plumedLabelWrong "\v<LABEL\=\@[^ #]*" contained
highlight link plumedLabelWrong Error

syntax region  plumedComment start="\v^\s*ENDPLUMED>" end="\%$" fold
syntax match   plumedComment excludenl "\v#.*$"
highlight link plumedComment Comment
endfunction

call PlumedDefineSyntax()


fun! PlumedGuessRegion()
" this is to find the match
" first, sync syntax
            syn sync fromstart
" find the syntactic attribute of the present region
            let col=col(".")
            let line=line(".")
            let key=""
            let stack=synstack(line,col)
            if(len(stack)>0)
              let key = synIDattr(stack[0], "name")
            endif
            if(key=~"^plumed[LC]*Line.*")
              return substitute(key,"^plumed[LC]*Line","","")
            endif
            return ""
endfunction

fun! PlumedContextManual()
  if(exists("b:plumed_helpfile"))
    quit
    return
  endif
  let m=PlumedGuessRegion()
  if(m=="")
    return
  else
    let name=s:path . "/help/" . m . ".txt"
    if(exists("b:plumed_helpfile_vertical"))
      execute 'rightbelow vsplit | view ' name
    else
      execute 'rightbelow split | view ' name
    endif
    let b:plumed_helpfile=1
" this is to allow closing the window with F2
    if exists("g:plumed_shortcuts")
      noremap  <buffer> <F2> :PHelp<CR>
    endif
  endif
endfunction

fun! PlumedManualV()
  let b:plumed_helpfile_vertical=1
endfunction

fun! PlumedManualH()
  unlet b:plumed_helpfile_vertical
endfunction

command! -nargs=0 PHelp call PlumedContextManual()

" autocomplete function
fun! PlumedComplete(findstart, base)
" this is to find the start of the word to be completed
          if a:findstart
            " locate the start of the word
            let line = getline('.')
            let start = col('.') - 1
            while start > 0 && line[start - 1] =~ '[a-zA-Z\_\=\-\.]'
              let start -= 1
            endwhile
            return start
          else
" this is to find the match
" first, sync syntax
            syn sync fromstart
" find the syntactic attribute of the present region
            let col=col(".")
            let line=line(".")
            let key=""
            if col!=1
              let key = synIDattr(synID(line, col-1, 1), "name")
            else
              let stack=synstack(line,col)
              if(len(stack)>0)
                let key = synIDattr(stack[0], "name")
              endif
            endif
            let comp=[]
" retrieve action name
" normalize ___ to -
            let key1=substitute(substitute(key,"^plumed[LC]*Line","",""),"___","-","g")
            if key ==""
" if outside of any region, complete with list of actions
              let comp=b:plumedActions
            elseif has_key(b:plumedDictionary,key1)
" if inside a region in the form "plumedLineXXX"
" complete with keywords associated to action XXX
              let comp=b:plumedDictionary[key1]
            endif
            " find months matching with "a:base"
            let res = []
            for m in comp
" this is to allow m to be a dictionary
" with a word and a one-liner
              if(type(m)==type({}))
                let n=m["word"]
              else
                let n=m
              endif
              if n =~ '^' . a:base
                if(n!="LABEL=" || key =~ "^plumedLine.*" || key =~ "^plumedCLine.*")
                call add(res, m)
                endif
              endif
" in principle comp could be a heterogeneous list
" so it should be unlet to iterate the loop
              unlet m
            endfor
"           if("..." =~ '^' . a:base && (key=~"^plumedLLine.*" || key=~"^plumedLine.*"))
"             call add(res,{"word":"...","menu":"(start multiline statement)"})
"           endif
"           if("..." =~ '^' . a:base && (key=~"^plumedLCLine.*" || key=~"^plumedCLine.*") && getline('.')=~'^\s*$')
"              call add(res,{"word":"...","menu":"(end multiline statement)"})
"           endif
            if("#" =~ '^' . a:base && key!="plumedComment") 
               call add(res,{"word":"#","menu":"(add comment)"})
            endif
            if("ENDPLUMED" =~ '^' . a:base && key =="")
               call add(res,{"word":"ENDPLUMED","menu":"(end input)"})
            endif
            return res
          endif
        endfun
setlocal omnifunc=PlumedComplete

" inspect the entire file to find lines containing
" non highlighted characters
fun! PlumedAnnotateSyntax()
" buffer where errors are written
  let buffer=[]
  let l=1
" loop over lines
  while l <= line("$")
    let line=getline(l)
    let p=0
" line is assumed right a priori
    let wrong=0
" position the cursor and redraw the screen
    call cursor(l,1)
    redraw! "! is required for some reason
    while p <len(line)
      let stack=synstack(l,p+1)
      if line[p] !~ "[ \t]"
        if(len(stack)==0)
          let wrong=1
        elseif(synIDattr(stack[len(stack)-1],"name")=~"^plumed[LC]*Line.*")
          let wrong=1
        endif
      endif
      let annotation=""
      for s in stack
        let annotation=annotation."+".synIDattr(s,"name")
      endfor
      call add(buffer,printf("ANNOTATION %5d %3d %s %s",l,p,line[p],annotation))
      let p=p+1
    endwhile
    
    if(wrong)
      call add(buffer,"ERROR AT LINE ".l." : ".line)
    endif
    let l=l+1
  endwhile
" dump the buffer on a new window
  new
  for l in buffer
    put=l
  endfor
endfun

