;+
; NAME:              BVLS_PS
;
; PURPOSE:           Bounded value least-square minimization (BVLS)
;
; USAGE:
;                    bvls_ps, a, b, bnd, x, istate, loopa, 
;                       RNORM=rnorm
;                       ITMAX=itmax,   
;                       EPS=eps       
;                       STATUS=status, /QUIET, /DOUBLE, 
;
; ARGUMENTS:
;       A   [input]
;          On entry A contains the M by N matrix, A.
;       B   [input]
;          On entry B contains the M-vector, B.
;       BND [input]
;          BND[0,J] is the lower bound for X[J].
;          BND[1,J] is the upper bound for X[J].
;          Require:  BND[0,J]  <=  BND[1,J].
;       X   [output]
;          Solution N-vector.
;       ISTATE [optional input/output] 
;          (N+1)-vector containing the state of each variable. If it is
;          given in input, it indicates a "warm" start, i.e. the program
;          will start from a guess of which variables are active (i.e.
;          strictly within their bounds). In output, it contains the final 
;          state.
;          This may be used when a series of related problems is treated
;          (pass the vector computed at a preious call for the next).
;          istate[n] is the total number of components at their bounds.  
;          The absolute values of the first nbound entries of  istate  are the 
;          indices of these `bound' components of  x. The sign of istate[j]
;          j=0...nbound-1, indicates whether  x[|istate[j]|-1] is at its upper 
;          or lower bound.  istate[j] is positive if the component is at its 
;          upper bound, negative if the component is at its lower bound.
;          istate[j], j=nbound,...,n-1  contain the indices of the active
;          components of x.
;
; KEYWORDS:   [Optional input]
;       LOOPA=loopa
;          number of iterations taken in the main loop, Loop A.
;       ITMAX=itmax,   
;          Maximum number of iterations permitted. Defaulted to 3*N, this is 
;          usually larger than required.  
;       EPS=eps,       
;          Determines the relative linear dependence of a column vector for a 
;          variable moved from its initial value.  The default value is
;          EPS=(machar()).eps (single precision).  
;       /DOUBLE,
;          Setting this keyword forces some computations to be performed in
;          double precision. This may be required if: (i) N is large (>10^6)
;          or (ii) if the A or B values span a large amplitude (e.g. 10^20).
;          This slows the computation, and is not required in most cases.
;       /QUIET
;          Verbosity control. By default the procedure prints the messages.
;
; KEYWORDS:   [Optional output]
;       RNORM=rnorm,  
;         Euclidean norm of the residual vector, || b - A*X ||.
;       STATUS=status,  
;         Indicates status on return.
;            0   Solution completed.
;            1   A is not a 2D array
;            2   B is not a vector
;            3   Sizes of A and B and not compatible
;            4   BND is not a 2D array
;            5   Size of BND is not consistent with B
;            6   Some lower bounds are excessively large
;            7   Some upper bounds are excessively small
;            8   Inconsistent bounds 
;            9   No free variable
;           10   Too many active variables in input STATE
;           11   Too many free variables
;           12   Failed to converge (reached ITMAX)
;
; DESCRIPTION:   
;   BVLS_PS solves linear least-squares problems with upper and lower bounds 
;   on the variables, using an active set strategy. It solves the problem:
;
;          min  || a.x - b ||     such that   bl <= x <= bu
;                            2
;     where
;               x  is an unknown n-vector
;               a  is a given m by n matrix
;               b  is a given  m-vector
;               bl is a given n-vector of lower bounds on the components of x.
;               bu is a given n-vector of upper bounds on the components of x.
;
;   The unconstrained least-squares problems for each candidate set of free 
;   variables are solved using the QR decomposition. BVLS has a "warm-start" 
;   feature permitting some of the variables to be initialized at their upper 
;   or lower bounds, which speeds the solution of a sequence of related 
;   problems.
;
;   See the article ``Bounded Variable Least Squares:  An Algorithm and
;   Applications'' by P.B. Stark and R.L. Parker, in the journal:
;   Computational Statistics, vol.10(2), (1995) for further description and
;   applications to minimum l-1, l-2 and l-infinity fitting problems, as well
;   as finding bounds on linear functionals subject to bounds on variables and
;   fitting linear data within l-1, l-2 or l-infinity measures of misfit.
;
; AUTHOR: 
;   Robert L. Parker and Philip B. Stark   (Version 3/19/90)
;   Copyright of this software is reserved by the authors; however, this
;   algorithm and subroutine may be used freely for non-commercial
;   purposes, and may be distributed freely for non-commercial purposes.
;   The authors do not warrant this software in any way: use it at your
;   own risk.
;
; HISTORY: 
;   The original Fortran77 code, can be found in the authors web
;   pages: http://www.stat.berkeley.edu/~stark/ 
;
;   It was converted to Fortran90 by Alan Miller (2000-11-23  Time: 23:41:42)
;   This can be found from: http://jblevins.org/mirror/amiller/
;
;   The GDL/IDL transcoding of this latter was completed by Philippe Prugniel 
;   on 2011 January 18th. It use LA_LEAST_SQUARES (IDL) to perform the QR 
;   factorization.
;   This routine is used in ULySS for replace the Lawson and Hanson version,
;   BVLS_LH. It is about 3 times faster, and using the warm start capability
;   provides another factor 3 gain when a minimization with several components
;   is made.
;-
; CATEGORY:    ULY_UTIL
;------------------------------------------------------------------------------

;--------------------Bounded Variable Least Squares---------------------
;        
;  Robert L. Parker                           Philip B. Stark
;  Scripps Institution of Oceanography        Department of Statistics
;  University of California, San Diego        University of California
;  La Jolla CA 92093                          Berkeley CA 94720-3860
;  rlparker@ucsd.edu                          stark@stat.berkeley.edu
;
;  Method: active variable method along the general plan of NNLS by
;  Lawson & Hanson, "Solving Least Squares Problems", Prentice-Hall 1974.
;  See Algorithm 23.10.  Step numbers in comment statements refer to their
;  scheme.
;  For more details and further uses, see the article
;  "Bounded Variable Least Squares:  An Algorithm and Applications"
;  by Stark and Parker in 1995 Computational Statistics.
;
;  A number of measures are taken to enhance numerical reliability:
;
; 1. As noted by Lawson and Hanson, roundoff errors in the computation of the
;   gradient of the misfit may cause a component on the bounds to appear to
;   want to become active, yet when the component is added to the active set,
;   it moves away from the feasible region.  In this case the component is not
;   made active, the gradient of the misfit with respect to a change in that
;   component is set to zero, and the program returns to the Kuhn-Tucker test.
;   Flag  ifrom5  is used in this test, which occurs at the end of Step 6.
;
;
; 2. When the least-squares minimizer after Step 6 is infeasible, it is used
;   in a convex interpolation with the previous solution to obtain a feasible
;   vector.  The constant in this interpolation is supposed to put at least
;   one component of  x  on a bound.  There can be difficulties:
;
; 2a. Sometimes, due to roundoff, no interpolated component ends up on
;   a bound.  The code in Step 11 uses the flag  jj, computed in Step 8,
;   to ensure that at least the component that determined the
;   interpolation constant  alpha  is moved to the appropriate bound.
;   This guarantees that what Lawson and Hanson call `Loop B' is finite.
;
; 2b. The code in Step 11 also incorporates Lawson and Hanson's feature that
;   any components remaining infeasible at this stage (which must be due to
;   roundoff) are moved to their nearer bound.
;
;
; 3. If the columns of A passed to QR are linearly dependent, the new
;   potentially active component is not introduced: the gradient of the
;   misfit with respect to that component is set to zero, and control
;   returns to the Kuhn-Tucker test.
;
;
; 4. When some of the columns of A are approximately linearly dependent,
;   we have observed cycling of active components: a component just moved
;   to a bound desires immediately to become active again; QR allows it
;   to become active and a different component is moved to its bound.
;   This component immediately wants to become active, which QR allows, and
;   the original component is moved back to its bound.  We have taken two
;   steps to avoid this problem:
;
; 4a. First, the column of the matrix  a  corresponding to the new
;   potentially active component is passed to QR as the last column of
;   its matrix.  This ordering tends to make a component recently moved
;   to a bound fail the test mentioned in (1), above.
;
; 4b. Second, we have incorporated a test that prohibits short cycles.
;   If the most recent successful change to the active set was to move
;   the component x(jj) to a bound, x(jj) is not permitted to reenter
;   the solution at this stage.  This test occurs just after checking
;   the Kuhn-Tucker conditions, and uses the flag  jj, set in Step 8.
;   The flag jj is reset after Step 6 if Step 6 was entered from
;   Step 5 indicating that a new component has successfully entered the
;   active set.  The test for resetting  jj  uses the flag  ifrom5,
;   which will not equal zero in case Step 6 was entered from Step 5.

pro bvls_ps, a, b, bnd, x, istate, $
             RNORM=rnorm,    $
             DOUBLE=double,  $
             EPS=eps,        $
             ITMAX=itmax,    $
             LOOP=loopa,     $
             STATUS=status,  $
             QUIET=quiet      

;----------------------First Executable Statement-----------------------
;  Step 1.  Check the dimension of the input arrays, initial values, etc.
;  Initialize everything--active and bound sets

status = 0
if n_elements(eps) eq 0 then eps = (machar()).eps   

;  Even for the double precision computations, avoid some unreasonable values
huge = (machar()).xmax             ;  single precision largest usable number

sz = size(a)
if sz[0] ne 2 then begin
   status = 1
   if not keyword_set(quiet) then message, /INFO, 'Matrix A should be a 2D array'
   return
endif
m = sz[1]
n = sz[2]
sz = size(b)
if sz[0] ne 1 then begin
   status = 2
   if not keyword_set(quiet) then message, /INFO, 'B should be a vector'
   return
endif
if sz[1] ne m then begin
   status = 3
   if not keyword_set(quiet) then $
      message, /INFO, 'The length of B ('+strtrim(sz[1],2)+ $
               ') is not consitent with A ('+strtrim(m,2)+')'
   return
endif
sz = size(bnd)
if sz[0] ne 2 then begin
   status = 4
   if not keyword_set(quiet) then message, /INFO, 'BND should be 2D array'
   return
endif
if sz[1] ne 2 then begin
   status = 5
   if not keyword_set(quiet) then message, /INFO, 'BND should be a 2xN array'
   return
endif
if sz[2] ne n then  begin
   status = 5
   if not keyword_set(quiet) then message, /INFO, 'BND should be a 2xN array'
   return
endif

if n_elements(itmax) eq 0 then itmax = 3 * n

bl = bnd[0,*]
nl = where(finite(bl) eq 0, cnt)  ; Take non-finite values as 'No bound'
if cnt gt 0 then bl[nl] = -huge
bu = bnd[1,*]
nl = where(finite(bu) eq 0, cnt)  ; Take non-finite values as 'No bound'
if cnt gt 0 then bu[nl] = huge

key =  n_elements(istate) eq n+1 ? 1 : 0

; If any of A and B is double, we should make the computation in double
if keyword_set(double) then type = size(1d, /TYP) $
else type = size(A, /TYP) > size(B, /TYP)

;  Initialize flags, etc.
mm = m < n
jj = -1
ifrom5 = 0
iact = 0

; Declaration of variables
x = make_array(n, TYPE=type)               ; solution to be returned
w = make_array(n, TYPE=type)

; The QR solver, la_least_squares, require A to be transposed
; If we favor less memory, we shall to the transpose when calling QR
; but for efficiency, we prefer to transpose A for once:
atr = transpose(a)

; No bl should be greater than huge nor bu smaller than -huge
j = where(bl ge huge, cnt)
if cnt gt 0 then begin
   status = 6
   if not keyword_set(quiet) then message, /INFO, 'Some lower bounds are excessively large'
   return
endif
j = where(bu le -huge, cnt)
if cnt gt 0 then begin
   status = 7
   if not keyword_set(quiet) then message, /INFO, 'Some upper bounds are excessively small'
   return
endif

;  Check consistency of given bounds  bl, bu.
j = where(bl gt -huge and bu lt huge, cnt)
if cnt gt 0 then begin
   bdiff = max(bu[j]-bl[j], MIN=bdifmi)
   if bdifmi lt 0 then begin
      status = 8
      if not keyword_set(quiet) then  message, /INFO, 'Inconsistent bounds'
      return
   endif
   if bdiff eq 0 then begin
      status = 9
      if not keyword_set(quiet) then  message, /INFO, 'No free variables .. check input bounds.'
      return
   endif
endif

;  In a fresh initialization (key = 0) bind all variables at their lower bounds.
;  If (key != 0), use the supplied  istate  vector to initialize the variables.
;  istate(n) contains the number of bound variables.  The absolute values of
;  the first nbound entries of  istate  are the indices of the bound variables.
;  The sign of each entry determines whether the indicated variable is at its 
;  upper (positive) or lower (negative) bound.
if key eq 0 then begin
  nbound = n
  nact = 0
  istate = -lindgen(nbound+1) - 1
  j = where(bl le -huge, cnt) 
  if cnt gt 0 then istate[j] *= -1    ; no lower bounds
  j = where(istate[0:n-1] gt 0 and bu ge huge, cnt, COMPLEM=isbound) 
  if cnt gt 0 then begin
     nbound -= cnt
     j = istate[j]
     if cnt lt n then istate[0:nbound-1] = istate[isbound]
     istate[nbound:n-1] = j           ; unconstrained
     istate[n] = nbound
  endif
endif else nbound = istate[n]
istate[n] = 0

nact = n - nbound
if nact gt mm then begin
   status = 10
   if not keyword_set(quiet) then  message, /INFO, 'Too many active variables in starting solution!'
   return
endif

; no loop but longer and not clearer (the perforance incidence is tiny)
;j1 = where(istate lt 0, cnt, COMPLEM=j2)
;if cnt gt 0 then begin
;   j = istate[j1]   
;   x[j] = bl[j]
;endif
;if cnt lt n then begin
;   j = istate[j2]
;   x[j] = bu[j]
;endif
for k=0, nbound-1 do begin
   j = ABS(istate[k]) - 1
   if istate[k] lt 0 then x[j] = bl[j] else x[j] = bu[j]
endfor

;  In a warm start (key != 0) initialize the active variables to (bl+bu)/2.
;   This is needed in case the initial QR results in active variables
;   out-of-bounds and Steps 8-11 get executed the first time through.
if nbound lt n then begin
   kk = istate[nbound:n-1] - 1
   j = where(bu[kk] lt huge and bl[kk] gt -huge, cnt, COMPLEM=k)
   if cnt gt 0 then x[kk[j]] = (bu[kk[j]]+bl[kk[j]]) / 2 $
   else if cnt lt n-nbound then begin
      kk = kk[k]
      j1 = where(bu[kk] lt huge, cnt2, COMPLEM=j2)
      if cnt2 gt 0 then x[kk[j1]] = bu[kk[j1]] $
      else if cnt2 lt cnt then begin
         kk = kk[j2]
         j3 = where(bl[kk] gt -huge, cnt3, COMPLEM=j4)
         if cnt3 gt 0 then x[kk[j3]] = bl[kk[j3]] $
         else if cnt3 lt cnt2 then x[kk[j4]] = 0
      endif
   endif
endif

key = nbound lt n ? 1 : 0          ; In cass some variables are unconstrained
;
;!  Compute bnorm, the norm of the data vector b, for reference.
bsq = total(b^2)
bnorm = sqrt(bsq)

;-----------------------------Main Loop---------------------------------
;  Initialization complete.  Begin major loop (Loop A).

for loopa=1,3*n do begin
;  Step 2.
;  Compute the residual vector b-a.x , the negative gradient vector
;   w(*), and the current objective value obj = || a.x - b ||.
;   The residual vector is stored in the mm+1'st column of act(*,*).

; no gain to make this following computation
;if loopa eq 1 then rsid = b - a # x  $ ; residuals
;else if nbound lt n then begin  
;    j = istate[nbound:n-1] - 1
;    rsid -= a[*,j] # x[j] 
;endif
;rsid1 = rsid ; for test

   rsid = b - a # x       ; residuals
;rsid2 = rsid ; for test
   obj = total(rsid^2)    
   w = a ## rsid          

  
;  Converged?  Stop if the misfit << || b ||, or if all components are
;   active (unless this is the first iteration from a warm start).
  IF SQRT(obj) le bnorm*eps OR (loopa gt 1 AND nbound eq 0) THEN begin
     istate[n] = nbound
     w[0] = SQRT(obj)
     rnorm = w[0]
     return
  endif
  
;  Add the contribution of the active components back into the residual.
  if nbound lt n then begin  
     j = istate[nbound:n-1] - 1
     rsid += a[*,j] # x[j] 
  endif
;rsid3 = rsid ; for test
  
;  The first iteration in a warm start requires immediate QR.
  IF loopa eq 1 AND key ne 0 then GOTO, jump150
;  
;  Steps 3, 4.
;  Find the bound element that most wants to be active.
  jump120: 

  worst = 0.
  if nbound gt 0 then begin
     ks = ABS(istate[0:nbound-1]) - 1
     bad = w[ks]
     neg = where(istate[0:nbound-1] lt 0, cnt)
     if cnt gt 0 then bad[neg] *= -1
     worst = min(bad, it)
     iact = ks[it]
  endif

;  Test whether the Kuhn-Tucker condition is met.
  if worst ge 0d then begin
     istate[n] = nbound
     w[0] = sqrt(obj)
     rnorm = w[0]
     return
  endif
  
;  the component  (iact)  is the one that most wants to become active.
;   If the last successful change in the active set was to move x(iact)
;   to a bound, don't let x(iact) in now: set the derivative of the misfit
;   with respect to x(iact) to zero and return to the Kuhn-Tucker test.
  IF iact eq jj THEN begin
    w[jj] = 0d
    GOTO, jump120
  ENDIF
  
;  Step 5.
;  Undo the effect of the new (potentially) active variable on the
;   residual vector.
  IF istate[it] gt 0 then bound = bu[iact]
  IF istate[it] lt 0 then bound = bl[iact]
  rsid += bound * a[*,iact]
;rsid4 = rsid ; for test
  
;  Set flag ifrom5, indicating that Step 6 was entered from Step 5.
;   This forms the basis of a test for instability: the gradient calculation
;   shows that x(iact) wants to join the active set; if QR puts x(iact) beyond
;   the bound from which it came, the gradient calculation was in error and
;   the variable should not have been introduced.
  ifrom5 = istate[it]
  
;  Swap the indices (in istate) of the new active variable and the
;   rightmost bound variable; `unbind' that location by decrementing nbound.
  istate[it] = istate[nbound-1]
  nbound -= 1
  nact += 1
  istate[nbound] = iact + 1
  
  if mm lt nact then begin
     status = 11
     if not keyword_set(quiet) then  message, /INFO, ' Too many free variables!'
     return
  endif
  
  jump150: 
;  Step 6.
;  Select the appropriate columns of  a  for QR.  For added stability, reverse 
;   the column ordering so that the most recent addition to the active set is 
;   in the last column.  


  if nact-n+nbound gt 0 then message, 'care ... not tested'
  if nact-n+nbound gt 0 then $
     jcl = [jcl[0:nact-n+nbound-1],reverse(istate[nbound:n-1]) - 1] $
  else jcl = reverse(istate[nbound:n-1]) - 1

  resq = m ge n ? 0 : -1

  stts = check_math(MASK=32)
  zz = la_least_squares(atr[jcl,*], rsid, DOUBLE=double)
  if stts eq 0 then stts = check_math(MASK=32) ; dismiss possible underflows
  
;  Test for linear dependence in QR, and for an instability that moves the
;   variable just introduced away from the feasible region (rather than into
;   the region or all the way through it).
;   In either case, remove the latest vector introduced from the active set
;   and adjust the residual vector accordingly.
;   Set the gradient component (w(iact)) to zero and return to the Kuhn-Tucker
;   test.
  if resq lt 0d OR $
     (ifrom5 gt 0 AND bu[iact] lt huge AND zz[nact-1] gt bu[iact]) OR  $
     (ifrom5 lt 0 AND bl[iact] gt -huge AND zz[nact-1] lt bl[iact]) THEN begin
;     message, /INFO, 'Warning, this feature has not been tested yet';  ...apparently OK
     nbound ++
     if  x[iact]-bu[iact] lt 0 then istate[nbound-1] *= -1
     nact --
     rsid -=  x[iact] * a[*,iact]
     ifrom5 = 0
     w[iact] = 0d
     goto, jump120
  endif
  
;  If Step 6 was entered from Step 5 and we are here, a new variable
;   has been successfully introduced into the active set; the last
;   variable that was fixed at a bound is again permitted to become active.
  IF ifrom5 ne 0 then jj = -1
  ifrom5 = 0
  
;   Step 7.  Check for strict feasibility of the new QR solution.
  j = istate[nbound:nbound+nact-1] - 1
  rzz = reverse(zz[0:nact-1])
  jmp = where((rzz lt bl[j] and  bl[j] gt -huge) or $
              (rzz gt bu[j] and  bu[j] lt huge), jump)
  if jump gt 0 then k1 = jmp[0] + 1 

  if jump eq 0 then begin ;  New iterate is feasible; back to the top.
     x[j] = rzz

  endif else begin

;  Steps 8, 9.
     alpha = 2d
     alf = alpha
     for k=k1,nact do begin
        j = istate[k+nbound-1] - 1
        IF zz[nact-k] gt bu[j] and bu[j] lt huge then alf = (bu[j]-x[j]) / ( zz[nact-k]-x[j])
        IF zz[nact-k] lt bl[j] and bl[j] gt -huge  then alf = (bl[j]-x[j]) / ( zz[nact-k]-x[j])
        IF alf lt alpha THEN begin
           alpha = alf
           jj = j
           sj = zz[nact-k]-bl[j] ge 0 ? 1d : -1d
        endif
     endfor
  
;  Step 10
     j = istate[nbound:nbound+nact-1] - 1
     x[j] += alpha * (reverse(zz[0:nact-1])-x[j])
  
;  Step 11.
;  Move the variable that determined alpha to the appropriate bound.
;   (jj is its index; sj is + if zz(jj)> bu(jj), - if zz(jj)<bl(jj) ).
;   If any other component of  x  is infeasible at this stage, it must
;   be due to roundoff.  Bind every infeasible component and every
;   component at a bound to the appropriate bound.  Correct the
;   residual vector for any variables moved to bounds.  Since at least
;   one variable is removed from the active set in this step, Loop B
;   (Steps 6-11) terminates after at most  nact  steps.
     noldb = nbound
     for k=1,nact do begin
        j = istate[k+noldb-1] - 1
        IF (bu[j] lt huge and x[j] gt bu[j]) OR $
           (j eq jj AND sj gt 0d) THEN begin ;  Move x(j) to its upper bound.
           x[j] = bu[j]
           istate[k+noldb-1] = istate[nbound]
           istate[nbound] = j+1
           nbound ++
           rsid -= bu[j] * a[*,j]
        endif else if (bl[j] gt -huge and x[j] le bl[j]) OR $
           (j eq jj AND sj lt 0d) THEN begin ;  Move x(j) to its lower bound.
           x[j] = bl[j]
           istate[k+noldb-1] = istate[nbound]
           istate[nbound] = -j-1
           nbound ++
           rsid -= bl[j] * a[*,j]
        endif
     endfor
     nact = n - nbound
  
;  If there are still active variables left repeat the QR; if not,
;    go back to the top.
     IF nact gt 0 then GOTO, jump150
  endelse ; jump=1
;  
endfor

status = 12
if not keyword_set(quiet) then  message, /INFO, 'Failed to converge! '

end ;  bvls_ps

;
PRO Test_BVLS
  ncases = 1000
;  ncases = 5
  ncols = 10
;  ncols = 3


; Generate artificial data satisfying:  Y = A.beta + noise
seed = 1
beta = randomu(seed, ncols)
beta = 4.0 * (beta - 0.5)
a = randomu(seed, ncases,ncols)
y = dblarr(ncases)
for cas=0l,ncases-1 do begin
   e = randomu(seed)
   y[cas] = total(a[cas,*] * beta); + 0.1*(e - 0.5)
endfor
huge = (machar()).xmax 
undefine, istate
bnd = dblarr(2,ncols) 
bnd[0,*] += -1; + (machar()).xmax ;- !VALUES.D_INFINITY
bnd[1,*] += !VALUES.D_INFINITY

;print, 'y=', y
;print, 'a=', a
;print, 'y/a=', y/a

help, a, y
bvls_ps, a, y, bnd, x, istate, LOOP=loopa
bvls_ps, a, y, bnd, x, istate, LOOP=loopa

print, ' Column   Original beta   Solution'
for j=0,ncols-1 do $
   print, FORM='(i5, 2f14.3)', j+1, beta[j], x[j]
print, FORM='(a, i4)', ' No. of iterations = ', loopa

return
END ; PROGRAM Test_BVLS
