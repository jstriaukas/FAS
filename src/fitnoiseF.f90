! ------------------------------------------------------------------
! fitnoiseF.f90: fits LASSO effective noise
! ------------------------------------------------------------------
SUBROUTINE fitnoiseF(qlev, nquant, esim, numboot, nobs, nvars, x, y, pf, dfmax, pmax, &
& nlam, flmin, ulam, eps, isd, intr, maxit, nalam, b0, beta, ibeta, nbeta, &
& alam, npass, jerr, quant, hatlam) !enhat

  IMPLICIT NONE
  ! -------- INPUT VARIABLES -------- !
  INTEGER :: numboot, nobs, nvars, dfmax, pmax, nlam, nalam, isd, intr
  INTEGER :: npass, jerr, maxit
  INTEGER :: ibeta(pmax), nbeta(nlam)
  INTEGER :: nquant ! --- number of quantiles
  DOUBLE PRECISION :: flmin, eps, maxlam
  DOUBLE PRECISION :: x(nobs, nvars), y(nobs), esim(nobs,numboot)
  DOUBLE PRECISION :: pf(nvars)
  DOUBLE PRECISION :: beta(pmax, nlam), b0(nlam)
  DOUBLE PRECISION :: ulam(nlam), alam(nlam)
  DOUBLE PRECISION :: qlev(nquant), quant(nlam, nquant) ! ---  quantiles
  DOUBLE PRECISION :: hatlam(nquant) ! --- \hat{lambda}
  !DOUBLE PRECISION :: enhat(nquant,numboot)  ! --- distribution of effective noise
  ! -------- LOCAL DECLARATIONS -------- !
  INTEGER :: j, l, nk, ierr
  INTEGER, DIMENSION(:), ALLOCATABLE :: ju
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xmean, xnorm, maj
  ! -------- ALLOCATE VARIABLES -------- !
  ALLOCATE(ju(1:nvars), STAT=ierr)
  jerr = jerr + ierr
  ALLOCATE(xmean(1:nvars), STAT=ierr)
  jerr = jerr + ierr
  ALLOCATE(xnorm(1:nvars), STAT=ierr)
  jerr = jerr + ierr
  ALLOCATE(maj(1:nvars), STAT=ierr)
  jerr = jerr + ierr
  IF (jerr /= 0) RETURN
  CALL chkvars(nobs, nvars, x, ju)
  IF (MAXVAL(pf) <= 0.0D0) THEN
    jerr = 10000
    RETURN
  END IF
  pf = MAX(0.0D0, pf)
  ! -------------------- STANDARDIZE & COMPUTE MAJ --------------------- !
  CALL standard(nobs, nvars, x, ju, isd, intr, xmean, xnorm, maj)
  ! -------------------- COMPUTE LAMBDA --------------------- !
  IF (ulam(1) .EQ. -1.0D0) THEN
    CALL maxlambda(nvars, nobs, x, y, pf, maxlam)
    ulam(1) = maxlam
    DO j = 2, nlam !
        ulam(j) = maxlam + (maxlam/nlam - maxlam) * (j - 1) / (nlam - 1)
        !--equidistant in log space
        !tmp = LOG(maxlam) + (LOG(maxlam*flmin) - LOG(maxlam)) * (j - 1) / (nlam - 1)
        !ulam(j) = EXP(tmp)
    END DO
  END IF
! -------------------- CALL lassofitpathF --------------------- !
CALL lassofitpathqF(qlev, nquant, esim, numboot, maj, nobs, nvars, x, y, ju, pf, dfmax, &
 & pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, ibeta, nbeta, alam, &
  & npass, jerr, intr, quant, hatlam) !enhat

  IF (jerr > 0) RETURN ! CHECK ERROR AFTER CALLING FUNCTION
  ! ----------- TRANSFORM BETA BACK TO THE ORIGINAL SCALE ----------- !
  DO l = 1, nalam
    nk = nbeta(l)
    IF (isd == 1) THEN
      DO j = 1, nk
        beta(j,l) = beta(j,l)/xnorm(ibeta(j))
      END DO
    END IF
    IF (intr == 1) THEN
        b0(l) = b0(l) - DOT_PRODUCT(beta(1:nk,l),xmean(ibeta(1:nk)))
    END IF
  END DO
  DEALLOCATE(ju,xmean,xnorm,maj)
  RETURN
END SUBROUTINE fitnoiseF


! --------------------------------- lassofitpathqF --------------------------------- !
SUBROUTINE lassofitpathqF(qlev, nquant, esim, numboot, maj, nobs, nvars, x, y, ju, pf, dfmax, &
& pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, m, nbeta, alam, &
& npass, jerr, intr, quant, hatlam)!enhat
    ! xlo, xhi, nmiss
    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER :: mnl, numboot, nobs, nvars, dfmax, pmax, nlam, maxit, nalam, npass, jerr, intr
    INTEGER :: ju(nvars), m(pmax), nbeta(nlam)
    INTEGER :: nquant ! --- number of quantiles
    DOUBLE PRECISION :: eps
    DOUBLE PRECISION :: x(nobs, nvars), y(nobs), maj(nvars), esim(nobs,numboot)
    DOUBLE PRECISION :: pf(nvars)
    DOUBLE PRECISION :: beta(pmax, nlam), b0(nlam)
    DOUBLE PRECISION :: ulam(nlam), alam(nlam)
    DOUBLE PRECISION :: qlev(nquant), quant(nlam,nquant) ! --- quantile levels and quantiles
    DOUBLE PRECISION :: hatlam(nquant) ! --- \hat{lambda}
    !DOUBLE PRECISION :: enhat(nquant,numboot) ! --- distribution of effective noise
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER,  PARAMETER :: mnlam = 6
    DOUBLE PRECISION :: tmp, d, dif, oldb, u, v, al, flmin,  nobsd, tmpb
    DOUBLE PRECISION,  DIMENSION(:),  ALLOCATABLE :: b, oldbeta, r
    DOUBLE PRECISION,  DIMENSION(:,:),  ALLOCATABLE :: rhat, xrhat, en
    INTEGER :: k, j, l, vrg, ctr, ierr, ni, me, pln, ib, jb, t
    INTEGER,  DIMENSION(:),  ALLOCATABLE :: mm
    ! -------- ALLOCATE VARIABLES -------- !
    ALLOCATE(b(0:nvars), STAT=jerr)
    ALLOCATE(oldbeta(0:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE(mm(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE(r(1:nobs), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE(rhat(1:nobs,1:numboot), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE(xrhat(1:nvars,1:numboot), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE(en(1:nlam,1:numboot), STAT=ierr)
    jerr = jerr + ierr
    IF (jerr /= 0) RETURN
    ! ---------------- INITIALIZATION ---------------- !
    r = y
    b = 0.0D0
    oldbeta = 0.0D0
    m = 0
    mm = 0
    npass = 0
    ni = npass
    mnl = MIN(mnlam,nlam)
    ju = 0
    nobsd = 0.0D0
    DO ib = 1, nobs
        nobsd = nobsd + 1.0D0
    END DO
    quant = 0.0D0
    hatlam = -1.0D0
    t = 0
    !enhat = 0.0D0
    ! ----------------- LAMBDA LOOP (OUTMOST LOOP) ------------------- !
    DO l = 1, nlam
        al = ulam(l)
        ctr = 0
        ! ------------------ OUTER LOOP -------------------- !
        DO

            IF (intr == 1) oldbeta(0) = b(0)
            IF (ni > 0) oldbeta(m(1:ni)) = b(m(1:ni))
            pln = 0
            ! ----------------- MIDDLE LOOP -------------------- !
            DO
                npass = npass + 1
                dif = 0.0D0
                pln = pln + 1
                DO k = 1, nvars
                    IF (pln == 1) THEN
                        oldb = b(k)
                        DO ! BEGIN PROXIMAL COORDINATE DESCENT
                            u = maj(k) * b(k) + DOT_PRODUCT(r,x(:,k))/nobs
                            v = ABS(u) - al * pf(k)
                            IF (v > 0.0D0) THEN
                                tmp = SIGN(v,u)/maj(k)
                            ELSE
                                tmp = 0.0D0
                            END IF
                            d = tmp - b(k)
                            IF (d**2 < eps) EXIT
                            b(k) = tmp
                            r = r - x(:,k) * d
                        END DO ! END PROXIMAL GRADIENT DESCENT
                        d = b(k) - oldb
                        IF (ABS(d) > 0.0D0) THEN
                            dif = dif + maj(k) * d**2
                            IF (mm(k) == 0) THEN
                                ni = ni + 1
                                IF (ni > pmax) EXIT
                                mm(k) = ni
                                m(ni) = k ! RECORD ACTIVE VARIABLES
                            END IF
                        END IF
                        IF (ABS(b(k))>0.0D0) ju(k) = 1
                    ELSE
                            IF (ju(k) == 1) THEN
                                oldb = b(k)
                                DO ! BEGIN PROXIMAL GRADIENT DESCENT
                                    !u = 0.0D0
                                    u = DOT_PRODUCT(r,x(:,k))
                                    u = maj(k) * b(k) + u/nobs
                                    v = ABS(u) - al * pf(k)
                                    IF (v > 0.0D0) THEN
                                        tmp = SIGN(v,u)/maj(k)
                                    ELSE
                                        tmp = 0.0D0
                                    END IF
                                    d = tmp - b(k)
                                    IF (d**2 < eps) EXIT
                                    b(k) = tmp
                                    r = r - x(:,k) * d
                                END DO ! END PROXIMAL GRADIENT DESCENT
                                d = b(k) - oldb
                                IF (ABS(d) > 0.0D0) THEN
                                        dif = MAX(dif, maj(k) * d**2)
                                END IF
                            END IF
                    END IF
                END DO
                IF (ni > pmax) EXIT
                IF (intr == 1) THEN
                    oldb = b(0)
                    DO ! BEGIN GRADIENT DESCENT
                        d = SUM(r)/nobs
                        IF (d**2 < eps) EXIT
                        b(0) = b(0) + d
                        r = r - d
                    END DO ! END GRADIENT DESCENT
                    d = b(0) - oldb
                    IF (ABS(d) > 0.0D0) dif = MAX(dif, d**2)
                END IF
                IF (dif < eps) EXIT
            END DO ! ----------> END MIDDLE LOOP
            IF (ni > pmax) EXIT
            ! -------------- FINAL CHECK ---------------- !
            vrg = 1
            IF (intr == 1) THEN
                IF ((b(0) - oldbeta(0))**2 >= eps) vrg = 0
            END IF
            DO j = 1, ni
                IF ((b(m(j)) - oldbeta(m(j)))**2 >= eps) THEN
                    vrg = 0
                    EXIT
                END IF
            END DO
            IF (vrg == 1) EXIT
            ctr = ctr + 1
            IF (ctr > maxit) THEN
                jerr = - l
                RETURN
            END IF
        END DO ! -------> END OUTER LOOP
        IF (numboot .NE. -1) THEN
            ! -------> COMPUTE EFFECTIVE NOISE
            DO ib = 1, numboot
                DO jb = 1, nobs ! nobs x numboot
                    rhat(jb,ib) = r(jb)*esim(jb,ib)
                END DO
            END DO
            xrhat = MATMUL(TRANSPOSE(x), rhat)/nobsd
            xrhat = 2*ABS(xrhat)
            DO ib = 1, numboot
                tmpb = MAXVAL(xrhat(:,ib))
                en(l, ib) = tmpb
            END DO
            CALL QUANTILE(en(l,:), qlev, nquant, numboot, quant(l,:))
            DO ib = 1, nquant
                IF (quant(l,ib) >= 2*ulam(l)) THEN
                    IF (hatlam(ib) .EQ. -1.0D0) THEN
                        IF (l .EQ. 1) THEN
                            hatlam(ib) = quant(1,ib)
                            !enhat(ib,:) = en(1,:)
                        ELSE
                            hatlam(ib) = quant(l,ib)
                            !enhat(ib,:) = en(l,:)
                        END IF
                        t = t + 1
                        IF (t .EQ. nquant) EXIT
                    END IF
                END IF
            END DO
        END IF
        ! ----------- FINAL UPDATE & SAVE RESULTS ------------ !
        IF (ni > pmax) THEN
            jerr = - 10000 - l
            EXIT
        END IF
        IF (ni > 0) beta(1:ni,l) = b(m(1:ni))
        nbeta(l) = ni
        b0(l) = b(0)
        alam(l) = al
        nalam = l
        IF (l < mnl) CYCLE
        IF (flmin >= 1.0D0) CYCLE
        me = COUNT(ABS(beta(1:ni,l)) > 0.0D0)
        IF (me > dfmax) EXIT
        IF (t .EQ. nquant) EXIT
    END DO ! -------> END LAMBDA LOOP
    DEALLOCATE(b,oldbeta,r,mm,rhat,xrhat,en)
    RETURN
END SUBROUTINE lassofitpathqF

SUBROUTINE QUANTILE(x, probs, np, n, qs)
    ! Calculate quantiles using type 7 (default in R)
    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER :: np, n
    DOUBLE PRECISION :: x(n), probs(np), qs(np), xs(n)
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER :: i, j, lo, hi
    DOUBLE PRECISION :: index, nd, lod, hid, h
    nd = 0.0D0
    DO i = 1, n
        nd = nd + 1.0D0
    END DO
    DO i = 1, np
        index = 1.0D0 + (nd - 1) * probs(i)
        lo = FLOOR(index)
        lod = 0.0D0
        DO j = 1, lo
            lod = lod + 1.0D0
        END DO
        hi = CEILING(index)
        hid = 0.0D0
        DO j = 1, hi
            hid = hid + 1.0D0
        END DO
        xs = x
        !CALL SORT(xs,n)
        CALL QUICKSORT(xs, 1, n) ! quicksort is much faster than simple sort, especially for large xs
        qs(i) = xs(lo)
        h = index - lod
        qs(i) = (1 - h) * qs(i) + h * xs(hi)
    END DO
END SUBROUTINE QUANTILE

RECURSIVE SUBROUTINE QUICKSORT(a, first, last)
    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    DOUBLE PRECISION ::  a(*), x, t
    INTEGER :: first, last
    INTEGER ::  i, j

    x = a( (first+last) / 2 )
    i = first
    j = last
    DO
    DO WHILE (a(i) < x)
        i=i+1
    END DO
    DO WHILE (x < a(j))
        j=j-1
    END DO
    IF (i >= j) EXIT
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
    END DO
    IF (first < i-1) CALL QUICKSORT(a, first, i-1)
    IF (j+1 < last)  CALL QUICKSORT(a, j+1, last)
END SUBROUTINE QUICKSORT

SUBROUTINE SORT(x, n)
    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER :: n
    DOUBLE PRECISION :: x(n)
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER :: K, L
    DOUBLE PRECISION :: T

    DO K = 1, n-1
        DO L = K+1, n
            IF (x(K) .GT. x(L)) THEN
                T = x(K)
                x(K) = x(L)
                x(L) = T
            END IF
        END DO
    END DO
END SUBROUTINE SORT
