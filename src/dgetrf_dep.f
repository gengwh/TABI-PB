*> \brief \b DGETRF
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DGETRF + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetrf.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetrf.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetrf.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGETRF computes an LU factorization of a general M-by-N matrix A
*> using partial pivoting with row interchanges.
*>
*> The factorization has the form
*>    A = P * L * U
*> where P is a permutation matrix, L is lower triangular with unit
*> diagonal elements (lower trapezoidal if m > n), and U is upper
*> triangular (upper trapezoidal if m < n).
*>
*> This is the right-looking Level 3 BLAS version of the algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the M-by-N matrix to be factored.
*>          On exit, the factors L and U from the factorization
*>          A = P*L*U; the unit diagonal elements of L are not stored.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (min(M,N))
*>          The pivot indices; for 1 <= i <= min(M,N), row i of the
*>          matrix was interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*>                has been completed, but the factor U is exactly
*>                singular, and division by zero will occur if it is used
*>                to solve a system of equations.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup doubleGEcomputational
*
*  =====================================================================
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGETRF2, DLASWP, DTRSM, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
         CALL DGETRF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        Use blocked code.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL DGETRF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of DGETRF
*
      END




*> \brief \b DGETRS
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DGETRS + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetrs.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetrs.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetrs.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          TRANS
*       INTEGER            INFO, LDA, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGETRS solves a system of linear equations
*>    A * X = B  or  A**T * X = B
*> with a general N-by-N matrix A using the LU factorization computed
*> by DGETRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          Specifies the form of the system of equations:
*>          = 'N':  A * X = B  (No transpose)
*>          = 'T':  A**T* X = B  (Transpose)
*>          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of columns
*>          of the matrix B.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          The factors L and U from the factorization A = P*L*U
*>          as computed by DGETRF.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*>          matrix was interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
*>          On entry, the right hand side matrix B.
*>          On exit, the solution matrix X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup doubleGEcomputational
*
*  =====================================================================
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASWP, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( NOTRAN ) THEN
*
*        Solve A * X = B.
*
*        Apply row interchanges to the right hand sides.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
*
*        Solve L*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve U*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
*
*        Solve A**T * X = B.
*
*        Solve U**T *X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve L**T *X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE,
     $               A, LDA, B, LDB )
*
*        Apply row interchanges to the solution vectors.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
*
      RETURN
*
*     End of DGETRS
*
      END



*> \brief \b DLASWP performs a series of row interchanges on a general rectangular matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DLASWP + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaswp.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaswp.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaswp.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*       .. Scalar Arguments ..
*       INTEGER            INCX, K1, K2, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLASWP performs a series of row interchanges on the matrix A.
*> One row interchange is initiated for each of rows K1 through K2 of A.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the matrix of column dimension N to which the row
*>          interchanges will be applied.
*>          On exit, the permuted matrix.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.
*> \endverbatim
*>
*> \param[in] K1
*> \verbatim
*>          K1 is INTEGER
*>          The first element of IPIV for which a row interchange will
*>          be done.
*> \endverbatim
*>
*> \param[in] K2
*> \verbatim
*>          K2 is INTEGER
*>          (K2-K1+1) is the number of elements of IPIV for which a row
*>          interchange will be done.
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (K1+(K2-K1)*abs(INCX))
*>          The vector of pivot indices. Only the elements in positions
*>          K1 through K1+(K2-K1)*abs(INCX) of IPIV are accessed.
*>          IPIV(K1+(K-K1)*abs(INCX)) = L implies rows K and L are to be
*>          interchanged.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          The increment between successive values of IPIV. If INCX
*>          is negative, the pivots are applied in reverse order.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup doubleOTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Modified by
*>   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      DOUBLE PRECISION   TEMP
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(K1+(I-K1)*abs(INCX)) for each of rows
*     K1 through K2.
*
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = K1 + ( K1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
*
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
*
      RETURN
*
*     End of DLASWP
*
      END




*> \brief \b DGETRF2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       RECURSIVE SUBROUTINE DGETRF2( M, N, A, LDA, IPIV, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGETRF2 computes an LU factorization of a general M-by-N matrix A
*> using partial pivoting with row interchanges.
*>
*> The factorization has the form
*>    A = P * L * U
*> where P is a permutation matrix, L is lower triangular with unit
*> diagonal elements (lower trapezoidal if m > n), and U is upper
*> triangular (upper trapezoidal if m < n).
*>
*> This is the recursive version of the algorithm. It divides
*> the matrix into four submatrices:
*>
*>        [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2
*>    A = [ -----|----- ]  with n1 = min(m,n)/2
*>        [  A21 | A22  ]       n2 = n-n1
*>
*>                                       [ A11 ]
*> The subroutine calls itself to factor [ --- ],
*>                                       [ A12 ]
*>                 [ A12 ]
*> do the swaps on [ --- ], solve A12, update A22,
*>                 [ A22 ]
*>
*> then calls itself to factor A22 and do the swaps on A21.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the M-by-N matrix to be factored.
*>          On exit, the factors L and U from the factorization
*>          A = P*L*U; the unit diagonal elements of L are not stored.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (min(M,N))
*>          The pivot indices; for 1 <= i <= min(M,N), row i of the
*>          matrix was interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*>                has been completed, but the factor U is exactly
*>                singular, and division by zero will occur if it is used
*>                to solve a system of equations.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup doubleGEcomputational
*
*  =====================================================================
      RECURSIVE SUBROUTINE DGETRF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   SFMIN, TEMP
      INTEGER            I, IINFO, N1, N2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      INTEGER            IDAMAX
      EXTERNAL           DLAMCH, IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DSCAL, DLASWP, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN

      IF ( M.EQ.1 ) THEN
*
*        Use unblocked code for one row case
*        Just need to handle IPIV and INFO
*
         IPIV( 1 ) = 1
         IF ( A(1,1).EQ.ZERO )
     $      INFO = 1
*
      ELSE IF( N.EQ.1 ) THEN
*
*        Use unblocked code for one column case
*
*
*        Compute machine safe minimum
*
         SFMIN = DLAMCH('S')
*
*        Find pivot and test for singularity
*
         I = IDAMAX( M, A( 1, 1 ), 1 )
         IPIV( 1 ) = I
         IF( A( I, 1 ).NE.ZERO ) THEN
*
*           Apply the interchange
*
            IF( I.NE.1 ) THEN
               TEMP = A( 1, 1 )
               A( 1, 1 ) = A( I, 1 )
               A( I, 1 ) = TEMP
            END IF
*
*           Compute elements 2:M of the column
*
            IF( ABS(A( 1, 1 )) .GE. SFMIN ) THEN
               CALL DSCAL( M-1, ONE / A( 1, 1 ), A( 2, 1 ), 1 )
            ELSE
               DO 10 I = 1, M-1
                  A( 1+I, 1 ) = A( 1+I, 1 ) / A( 1, 1 )
   10          CONTINUE
            END IF
*
         ELSE
            INFO = 1
         END IF
*
      ELSE
*
*        Use recursive code
*
         N1 = MIN( M, N ) / 2
         N2 = N-N1
*
*               [ A11 ]
*        Factor [ --- ]
*               [ A21 ]
*
         CALL DGETRF2( M, N1, A, LDA, IPIV, IINFO )

         IF ( INFO.EQ.0 .AND. IINFO.GT.0 )
     $      INFO = IINFO
*
*                              [ A12 ]
*        Apply interchanges to [ --- ]
*                              [ A22 ]
*
         CALL DLASWP( N2, A( 1, N1+1 ), LDA, 1, N1, IPIV, 1 )
*
*        Solve A12
*
         CALL DTRSM( 'L', 'L', 'N', 'U', N1, N2, ONE, A, LDA,
     $               A( 1, N1+1 ), LDA )
*
*        Update A22
*
         CALL DGEMM( 'N', 'N', M-N1, N2, N1, -ONE, A( N1+1, 1 ), LDA,
     $               A( 1, N1+1 ), LDA, ONE, A( N1+1, N1+1 ), LDA )
*
*        Factor A22
*
         CALL DGETRF2( M-N1, N2, A( N1+1, N1+1 ), LDA, IPIV( N1+1 ),
     $                 IINFO )
*
*        Adjust INFO and the pivot indices
*
         IF ( INFO.EQ.0 .AND. IINFO.GT.0 )
     $      INFO = IINFO + N1
         DO 20 I = N1+1, MIN( M, N )
            IPIV( I ) = IPIV( I ) + N1
   20    CONTINUE
*
*        Apply interchanges to A21
*
         CALL DLASWP( N1, A( 1, 1 ), LDA, N1+1, MIN( M, N), IPIV, 1 )
*
      END IF
      RETURN
*
*     End of DGETRF2
*
      END



*> \brief \b DGEMM
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA,BETA
*       INTEGER K,LDA,LDB,LDC,M,N
*       CHARACTER TRANSA,TRANSB
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGEMM  performs one of the matrix-matrix operations
*>
*>    C := alpha*op( A )*op( B ) + beta*C,
*>
*> where  op( X ) is one of
*>
*>    op( X ) = X   or   op( X ) = X**T,
*>
*> alpha and beta are scalars, and A, B and C are matrices, with op( A )
*> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANSA
*> \verbatim
*>          TRANSA is CHARACTER*1
*>           On entry, TRANSA specifies the form of op( A ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSA = 'N' or 'n',  op( A ) = A.
*>
*>              TRANSA = 'T' or 't',  op( A ) = A**T.
*>
*>              TRANSA = 'C' or 'c',  op( A ) = A**T.
*> \endverbatim
*>
*> \param[in] TRANSB
*> \verbatim
*>          TRANSB is CHARACTER*1
*>           On entry, TRANSB specifies the form of op( B ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSB = 'N' or 'n',  op( B ) = B.
*>
*>              TRANSB = 'T' or 't',  op( B ) = B**T.
*>
*>              TRANSB = 'C' or 'c',  op( B ) = B**T.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry,  M  specifies  the number  of rows  of the  matrix
*>           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry,  N  specifies the number  of columns of the matrix
*>           op( B ) and the number of columns of the matrix C. N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>           On entry,  K  specifies  the number of columns of the matrix
*>           op( A ) and the number of rows of the matrix op( B ). K must
*>           be at least  zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension ( LDA, ka ), where ka is
*>           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*>           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*>           part of the array  A  must contain the matrix  A,  otherwise
*>           the leading  k by m  part of the array  A  must contain  the
*>           matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*>           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*>           least  max( 1, k ).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension ( LDB, kb ), where kb is
*>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*>           part of the array  B  must contain the matrix  B,  otherwise
*>           the leading  n by k  part of the array  B  must contain  the
*>           matrix B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*>           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*>           least  max( 1, n ).
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION.
*>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*>           supplied as zero then C need not be set on input.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array, dimension ( LDC, N )
*>           Before entry, the leading  m by n  part of the array  C must
*>           contain the matrix  C,  except when  beta  is zero, in which
*>           case C need not be set on entry.
*>           On exit, the array  C  is overwritten by the  m by n  matrix
*>           ( alpha*op( A )*op( B ) + beta*C ).
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>           On entry, LDC specifies the first dimension of C as declared
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least
*>           max( 1, m ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup double_blas_level3
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 3 Blas routine.
*>
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC max
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NROWA,NROWB
      LOGICAL NOTA,NOTB
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
*     ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA and NROWB  as the number of rows of  A
*     and  B  respectively.
*
      nota = lsame(transa,'N')
      notb = lsame(transb,'N')
      IF (nota) THEN
          nrowa = m
      ELSE
          nrowa = k
      END IF
      IF (notb) THEN
          nrowb = k
      ELSE
          nrowb = n
      END IF
*
*     Test the input parameters.
*
      info = 0
      IF ((.NOT.nota) .AND. (.NOT.lsame(transa,'C')) .AND.
     +    (.NOT.lsame(transa,'T'))) THEN
          info = 1
      ELSE IF ((.NOT.notb) .AND. (.NOT.lsame(transb,'C')) .AND.
     +         (.NOT.lsame(transb,'T'))) THEN
          info = 2
      ELSE IF (m.LT.0) THEN
          info = 3
      ELSE IF (n.LT.0) THEN
          info = 4
      ELSE IF (k.LT.0) THEN
          info = 5
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 8
      ELSE IF (ldb.LT.max(1,nrowb)) THEN
          info = 10
      ELSE IF (ldc.LT.max(1,m)) THEN
          info = 13
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DGEMM ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((m.EQ.0) .OR. (n.EQ.0) .OR.
     +    (((alpha.EQ.zero).OR. (k.EQ.0)).AND. (beta.EQ.one))) RETURN
*
*     And if  alpha.eq.zero.
*
      IF (alpha.EQ.zero) THEN
          IF (beta.EQ.zero) THEN
              DO 20 j = 1,n
                  DO 10 i = 1,m
                      c(i,j) = zero
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 j = 1,n
                  DO 30 i = 1,m
                      c(i,j) = beta*c(i,j)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (notb) THEN
          IF (nota) THEN
*
*           Form  C := alpha*A*B + beta*C.
*
              DO 90 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 50 i = 1,m
                          c(i,j) = zero
   50                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 60 i = 1,m
                          c(i,j) = beta*c(i,j)
   60                 CONTINUE
                  END IF
                  DO 80 l = 1,k
                      temp = alpha*b(l,j)
                      DO 70 i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
   70                 CONTINUE
   80             CONTINUE
   90         CONTINUE
          ELSE
*
*           Form  C := alpha*A**T*B + beta*C
*
              DO 120 j = 1,n
                  DO 110 i = 1,m
                      temp = zero
                      DO 100 l = 1,k
                          temp = temp + a(l,i)*b(l,j)
  100                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  110             CONTINUE
  120         CONTINUE
          END IF
      ELSE
          IF (nota) THEN
*
*           Form  C := alpha*A*B**T + beta*C
*
              DO 170 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 130 i = 1,m
                          c(i,j) = zero
  130                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 140 i = 1,m
                          c(i,j) = beta*c(i,j)
  140                 CONTINUE
                  END IF
                  DO 160 l = 1,k
                      temp = alpha*b(j,l)
                      DO 150 i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
  150                 CONTINUE
  160             CONTINUE
  170         CONTINUE
          ELSE
*
*           Form  C := alpha*A**T*B**T + beta*C
*
              DO 200 j = 1,n
                  DO 190 i = 1,m
                      temp = zero
                      DO 180 l = 1,k
                          temp = temp + a(l,i)*b(j,l)
  180                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  190             CONTINUE
  200         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DGEMM
*
      END

*> \brief \b XERBLA
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE XERBLA( SRNAME, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER*(*)      SRNAME
*       INTEGER            INFO
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> XERBLA  is an error handler for the LAPACK routines.
*> It is called by an LAPACK routine if an input parameter has an
*> invalid value.  A message is printed and execution stops.
*>
*> Installers may consider modifying the STOP statement in order to
*> call system-specific exception-handling facilities.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SRNAME
*> \verbatim
*>          SRNAME is CHARACTER*(*)
*>          The name of the routine which called XERBLA.
*> \endverbatim
*>
*> \param[in] INFO
*> \verbatim
*>          INFO is INTEGER
*>          The position of the invalid parameter in the parameter list
*>          of the calling routine.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup aux_blas
*
*  =====================================================================
      SUBROUTINE xerbla( SRNAME, INFO )
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER*(*)      SRNAME
      INTEGER            INFO
*     ..
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          len_trim
*     ..
*     .. Executable Statements ..
*
      WRITE( *, fmt = 9999 )srname( 1:len_trim( srname ) ), info
*
      stop
*
 9999 FORMAT( ' ** On entry to ', a, ' parameter number ', i2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END




*> \brief \b DTRSM
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA
*       INTEGER LDA,LDB,M,N
*       CHARACTER DIAG,SIDE,TRANSA,UPLO
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),B(LDB,*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DTRSM  solves one of the matrix equations
*>
*>    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*>
*> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*>
*>    op( A ) = A   or   op( A ) = A**T.
*>
*> The matrix X is overwritten on B.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>           On entry, SIDE specifies whether op( A ) appears on the left
*>           or right of X as follows:
*>
*>              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*>
*>              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the matrix A is an upper or
*>           lower triangular matrix as follows:
*>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*> \endverbatim
*>
*> \param[in] TRANSA
*> \verbatim
*>          TRANSA is CHARACTER*1
*>           On entry, TRANSA specifies the form of op( A ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSA = 'N' or 'n'   op( A ) = A.
*>
*>              TRANSA = 'T' or 't'   op( A ) = A**T.
*>
*>              TRANSA = 'C' or 'c'   op( A ) = A**T.
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>           On entry, DIAG specifies whether or not A is unit triangular
*>           as follows:
*>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry, M specifies the number of rows of B. M must be at
*>           least zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of B.  N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*>           zero then  A is not referenced and  B need not be set before
*>           entry.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension ( LDA, k ),
*>           where k is m when SIDE = 'L' or 'l'
*>             and k is n when SIDE = 'R' or 'r'.
*>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*>           upper triangular part of the array  A must contain the upper
*>           triangular matrix  and the strictly lower triangular part of
*>           A is not referenced.
*>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*>           lower triangular part of the array  A must contain the lower
*>           triangular matrix  and the strictly upper triangular part of
*>           A is not referenced.
*>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*>           A  are not referenced either,  but are assumed to be  unity.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*>           then LDA must be at least max( 1, n ).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension ( LDB, N )
*>           Before entry,  the leading  m by n part of the array  B must
*>           contain  the  right-hand  side  matrix  B,  and  on exit  is
*>           overwritten by the solution matrix  X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least
*>           max( 1, m ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup double_blas_level3
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 3 Blas routine.
*>
*>
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE dtrsm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*
*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC max
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
*     ..
*
*     Test the input parameters.
*
      lside = lsame(side,'L')
      IF (lside) THEN
          nrowa = m
      ELSE
          nrowa = n
      END IF
      nounit = lsame(diag,'N')
      upper = lsame(uplo,'U')
*
      info = 0
      IF ((.NOT.lside) .AND. (.NOT.lsame(side,'R'))) THEN
          info = 1
      ELSE IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
          info = 2
      ELSE IF ((.NOT.lsame(transa,'N')) .AND.
     +         (.NOT.lsame(transa,'T')) .AND.
     +         (.NOT.lsame(transa,'C'))) THEN
          info = 3
      ELSE IF ((.NOT.lsame(diag,'U')) .AND. (.NOT.lsame(diag,'N'))) THEN
          info = 4
      ELSE IF (m.LT.0) THEN
          info = 5
      ELSE IF (n.LT.0) THEN
          info = 6
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 9
      ELSE IF (ldb.LT.max(1,m)) THEN
          info = 11
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DTRSM ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (m.EQ.0 .OR. n.EQ.0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (alpha.EQ.zero) THEN
          DO 20 j = 1,n
              DO 10 i = 1,m
                  b(i,j) = zero
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
*
*     Start the operations.
*
      IF (lside) THEN
          IF (lsame(transa,'N')) THEN
*
*           Form  B := alpha*inv( A )*B.
*
              IF (upper) THEN
                  DO 60 j = 1,n
                      IF (alpha.NE.one) THEN
                          DO 30 i = 1,m
                              b(i,j) = alpha*b(i,j)
   30                     CONTINUE
                      END IF
                      DO 50 k = m,1,-1
                          IF (b(k,j).NE.zero) THEN
                              IF (nounit) b(k,j) = b(k,j)/a(k,k)
                              DO 40 i = 1,k - 1
                                  b(i,j) = b(i,j) - b(k,j)*a(i,k)
   40                         CONTINUE
                          END IF
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 100 j = 1,n
                      IF (alpha.NE.one) THEN
                          DO 70 i = 1,m
                              b(i,j) = alpha*b(i,j)
   70                     CONTINUE
                      END IF
                      DO 90 k = 1,m
                          IF (b(k,j).NE.zero) THEN
                              IF (nounit) b(k,j) = b(k,j)/a(k,k)
                              DO 80 i = k + 1,m
                                  b(i,j) = b(i,j) - b(k,j)*a(i,k)
   80                         CONTINUE
                          END IF
   90                 CONTINUE
  100             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*inv( A**T )*B.
*
              IF (upper) THEN
                  DO 130 j = 1,n
                      DO 120 i = 1,m
                          temp = alpha*b(i,j)
                          DO 110 k = 1,i - 1
                              temp = temp - a(k,i)*b(k,j)
  110                     CONTINUE
                          IF (nounit) temp = temp/a(i,i)
                          b(i,j) = temp
  120                 CONTINUE
  130             CONTINUE
              ELSE
                  DO 160 j = 1,n
                      DO 150 i = m,1,-1
                          temp = alpha*b(i,j)
                          DO 140 k = i + 1,m
                              temp = temp - a(k,i)*b(k,j)
  140                     CONTINUE
                          IF (nounit) temp = temp/a(i,i)
                          b(i,j) = temp
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (lsame(transa,'N')) THEN
*
*           Form  B := alpha*B*inv( A ).
*
              IF (upper) THEN
                  DO 210 j = 1,n
                      IF (alpha.NE.one) THEN
                          DO 170 i = 1,m
                              b(i,j) = alpha*b(i,j)
  170                     CONTINUE
                      END IF
                      DO 190 k = 1,j - 1
                          IF (a(k,j).NE.zero) THEN
                              DO 180 i = 1,m
                                  b(i,j) = b(i,j) - a(k,j)*b(i,k)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
                      IF (nounit) THEN
                          temp = one/a(j,j)
                          DO 200 i = 1,m
                              b(i,j) = temp*b(i,j)
  200                     CONTINUE
                      END IF
  210             CONTINUE
              ELSE
                  DO 260 j = n,1,-1
                      IF (alpha.NE.one) THEN
                          DO 220 i = 1,m
                              b(i,j) = alpha*b(i,j)
  220                     CONTINUE
                      END IF
                      DO 240 k = j + 1,n
                          IF (a(k,j).NE.zero) THEN
                              DO 230 i = 1,m
                                  b(i,j) = b(i,j) - a(k,j)*b(i,k)
  230                         CONTINUE
                          END IF
  240                 CONTINUE
                      IF (nounit) THEN
                          temp = one/a(j,j)
                          DO 250 i = 1,m
                              b(i,j) = temp*b(i,j)
  250                     CONTINUE
                      END IF
  260             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*B*inv( A**T ).
*
              IF (upper) THEN
                  DO 310 k = n,1,-1
                      IF (nounit) THEN
                          temp = one/a(k,k)
                          DO 270 i = 1,m
                              b(i,k) = temp*b(i,k)
  270                     CONTINUE
                      END IF
                      DO 290 j = 1,k - 1
                          IF (a(j,k).NE.zero) THEN
                              temp = a(j,k)
                              DO 280 i = 1,m
                                  b(i,j) = b(i,j) - temp*b(i,k)
  280                         CONTINUE
                          END IF
  290                 CONTINUE
                      IF (alpha.NE.one) THEN
                          DO 300 i = 1,m
                              b(i,k) = alpha*b(i,k)
  300                     CONTINUE
                      END IF
  310             CONTINUE
              ELSE
                  DO 360 k = 1,n
                      IF (nounit) THEN
                          temp = one/a(k,k)
                          DO 320 i = 1,m
                              b(i,k) = temp*b(i,k)
  320                     CONTINUE
                      END IF
                      DO 340 j = k + 1,n
                          IF (a(j,k).NE.zero) THEN
                              temp = a(j,k)
                              DO 330 i = 1,m
                                  b(i,j) = b(i,j) - temp*b(i,k)
  330                         CONTINUE
                          END IF
  340                 CONTINUE
                      IF (alpha.NE.one) THEN
                          DO 350 i = 1,m
                              b(i,k) = alpha*b(i,k)
  350                     CONTINUE
                      END IF
  360             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of DTRSM
*
      END



*> \brief \b DLAMCH
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
*
*     .. Scalar Arguments ..
*     CHARACTER          CMACH
*     ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAMCH determines double precision machine parameters.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] CMACH
*> \verbatim
*>          CMACH is CHARACTER*1
*>          Specifies the value to be returned by DLAMCH:
*>          = 'E' or 'e',   DLAMCH := eps
*>          = 'S' or 's ,   DLAMCH := sfmin
*>          = 'B' or 'b',   DLAMCH := base
*>          = 'P' or 'p',   DLAMCH := eps*base
*>          = 'N' or 'n',   DLAMCH := t
*>          = 'R' or 'r',   DLAMCH := rnd
*>          = 'M' or 'm',   DLAMCH := emin
*>          = 'U' or 'u',   DLAMCH := rmin
*>          = 'L' or 'l',   DLAMCH := emax
*>          = 'O' or 'o',   DLAMCH := rmax
*>          where
*>          eps   = relative machine precision
*>          sfmin = safe minimum, such that 1/sfmin does not overflow
*>          base  = base of the machine
*>          prec  = eps*base
*>          t     = number of (base) digits in the mantissa
*>          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*>          emin  = minimum exponent before (gradual) underflow
*>          rmin  = underflow threshold - base**(emin-1)
*>          emax  = largest exponent before overflow
*>          rmax  = overflow threshold  - (base**emax)*(1-eps)
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
 
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      DOUBLE PRECISION FUNCTION dlamch( CMACH )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          cmach
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   one, zero
      parameter( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   rnd, eps, sfmin, small, rmach
*     ..
*     .. External Functions ..
      LOGICAL            lsame
      EXTERNAL           lsame
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          digits, epsilon, huge, maxexponent,
     $                   minexponent, radix, tiny
*     ..
*     .. Executable Statements ..
*
*
*     Assume rounding, not chopping. Always.
*
      rnd = one
*
      IF( one.EQ.rnd ) THEN
         eps = epsilon(zero) * 0.5
      ELSE
         eps = epsilon(zero)
      END IF
*
      IF( lsame( cmach, 'E' ) ) THEN
         rmach = eps
      ELSE IF( lsame( cmach, 'S' ) ) THEN
         sfmin = tiny(zero)
         small = one / huge(zero)
         IF( small.GE.sfmin ) THEN
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            sfmin = small*( one+eps )
         END IF
         rmach = sfmin
      ELSE IF( lsame( cmach, 'B' ) ) THEN
         rmach = radix(zero)
      ELSE IF( lsame( cmach, 'P' ) ) THEN
         rmach = eps * radix(zero)
      ELSE IF( lsame( cmach, 'N' ) ) THEN
         rmach = digits(zero)
      ELSE IF( lsame( cmach, 'R' ) ) THEN
         rmach = rnd
      ELSE IF( lsame( cmach, 'M' ) ) THEN
         rmach = minexponent(zero)
      ELSE IF( lsame( cmach, 'U' ) ) THEN
         rmach = tiny(zero)
      ELSE IF( lsame( cmach, 'L' ) ) THEN
         rmach = maxexponent(zero)
      ELSE IF( lsame( cmach, 'O' ) ) THEN
         rmach = huge(zero)
      ELSE
         rmach = zero
      END IF
*
      dlamch = rmach
      RETURN
*
*     End of DLAMCH
*
      END
************************************************************************
*> \brief \b DLAMC3
*> \details
*> \b Purpose:
*> \verbatim
*> DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*> the addition of  A  and  B ,  for use in situations where optimizers
*> might hold one of these in a register.
*> \endverbatim
*> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
*> \param[in] A
*> \verbatim
*>          A is a DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is a DOUBLE PRECISION
*>          The values A and B.
*> \endverbatim
*>
      DOUBLE PRECISION FUNCTION dlamc3( A, B )
*
*  -- LAPACK auxiliary routine --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   a, b
*     ..
* =====================================================================
*
*     .. Executable Statements ..
*
      dlamc3 = a + b
*
      RETURN
*
*     End of DLAMC3
*
      END
*
************************************************************************

*> \brief \b LSAME
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       LOGICAL FUNCTION LSAME(CA,CB)
*
*       .. Scalar Arguments ..
*       CHARACTER CA,CB
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> LSAME returns .TRUE. if CA is the same letter as CB regardless of
*> case.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] CA
*> \verbatim
*>          CA is CHARACTER*1
*> \endverbatim
*>
*> \param[in] CB
*> \verbatim
*>          CB is CHARACTER*1
*>          CA and CB specify the single characters to be compared.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup aux_blas
*
*  =====================================================================
      LOGICAL FUNCTION lsame(CA,CB)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER ca,cb
*     ..
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC ichar
*     ..
*     .. Local Scalars ..
      INTEGER inta,intb,zcode
*     ..
*
*     Test if the characters are equal
*
      lsame = ca .EQ. cb
      IF (lsame) RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      zcode = ichar('Z')
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      inta = ichar(ca)
      intb = ichar(cb)
*
      IF (zcode.EQ.90 .OR. zcode.EQ.122) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
          IF (inta.GE.97 .AND. inta.LE.122) inta = inta - 32
          IF (intb.GE.97 .AND. intb.LE.122) intb = intb - 32
*
      ELSE IF (zcode.EQ.233 .OR. zcode.EQ.169) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
          IF (inta.GE.129 .AND. inta.LE.137 .OR.
     +        inta.GE.145 .AND. inta.LE.153 .OR.
     +        inta.GE.162 .AND. inta.LE.169) inta = inta + 64
          IF (intb.GE.129 .AND. intb.LE.137 .OR.
     +        intb.GE.145 .AND. intb.LE.153 .OR.
     +        intb.GE.162 .AND. intb.LE.169) intb = intb + 64
*
      ELSE IF (zcode.EQ.218 .OR. zcode.EQ.250) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
          IF (inta.GE.225 .AND. inta.LE.250) inta = inta - 32
          IF (intb.GE.225 .AND. intb.LE.250) intb = intb - 32
      END IF
      lsame = inta .EQ. intb
*
*     RETURN
*
*     End of LSAME
*
      END



! *> \brief \b DSCAL
! *
! *  =========== DOCUMENTATION ===========
! *
! * Online html documentation available at
! *            http://www.netlib.org/lapack/explore-html/
! *
! *  Definition:
! *  ===========
! *
! *       SUBROUTINE DSCAL(N,DA,DX,INCX)
! *
! *       .. Scalar Arguments ..
! *       DOUBLE PRECISION DA
! *       INTEGER INCX,N
! *       ..
! *       .. Array Arguments ..
! *       DOUBLE PRECISION DX(*)
! *       ..
! *
! *
! *> \par Purpose:
! *  =============
! *>
! *> \verbatim
! *>
! *>    DSCAL scales a vector by a constant.
! *>    uses unrolled loops for increment equal to 1.
! *> \endverbatim
! *
! *  Arguments:
! *  ==========
! *
! *> \param[in] N
! *> \verbatim
! *>          N is INTEGER
! *>         number of elements in input vector(s)
! *> \endverbatim
! *>
! *> \param[in] DA
! *> \verbatim
! *>          DA is DOUBLE PRECISION
! *>           On entry, DA specifies the scalar alpha.
! *> \endverbatim
! *>
! *> \param[in,out] DX
! *> \verbatim
! *>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
! *> \endverbatim
! *>
! *> \param[in] INCX
! *> \verbatim
! *>          INCX is INTEGER
! *>         storage spacing between elements of DX
! *> \endverbatim
! *
! *  Authors:
! *  ========
! *
! *> \author Univ. of Tennessee
! *> \author Univ. of California Berkeley
! *> \author Univ. of Colorado Denver
! *> \author NAG Ltd.
! *
! *> \ingroup double_blas_level1
! *
! *> \par Further Details:
! *  =====================
! *>
! *> \verbatim
! *>
! *>     jack dongarra, linpack, 3/11/78.
! *>     modified 3/93 to return if incx .le. 0.
! *>     modified 12/3/93, array(1) declarations changed to array(*)
! *> \endverbatim
! *>
! *  =====================================================================
!       SUBROUTINE dscal(N,DA,DX,INCX)
! *
! *  -- Reference BLAS level1 routine --
! *  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
! *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
! *
! *     .. Scalar Arguments ..
!       DOUBLE PRECISION DA
!       INTEGER INCX,N
! *     ..
! *     .. Array Arguments ..
!       DOUBLE PRECISION DX(*)
! *     ..
! *
! *  =====================================================================
! *
! *     .. Local Scalars ..
!       INTEGER I,M,MP1,NINCX
! *     .. Parameters ..
!       DOUBLE PRECISION ONE
!       parameter(one=1.0d+0)
! *     ..
! *     .. Intrinsic Functions ..
!       INTRINSIC mod
! *     ..
!       IF (n.LE.0 .OR. incx.LE.0 .OR. da.EQ.one) RETURN
!       IF (incx.EQ.1) THEN
! *
! *        code for increment equal to 1
! *
! *
! *        clean-up loop
! *
!          m = mod(n,5)
!          IF (m.NE.0) THEN
!             DO i = 1,m
!                dx(i) = da*dx(i)
!             END DO
!             IF (n.LT.5) RETURN
!          END IF
!          mp1 = m + 1
!          DO i = mp1,n,5
!             dx(i) = da*dx(i)
!             dx(i+1) = da*dx(i+1)
!             dx(i+2) = da*dx(i+2)
!             dx(i+3) = da*dx(i+3)
!             dx(i+4) = da*dx(i+4)
!          END DO
!       ELSE
! *
! *        code for increment not equal to 1
! *
!          nincx = n*incx
!          DO i = 1,nincx,incx
!             dx(i) = da*dx(i)
!          END DO
!       END IF
!       RETURN
! *
! *     End of DSCAL
! *
!       END


*> \brief \b IDAMAX
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       INTEGER FUNCTION IDAMAX(N,DX,INCX)
*
*       .. Scalar Arguments ..
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION DX(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    IDAMAX finds the index of the first element having maximum absolute value.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         number of elements in input vector(s)
*> \endverbatim
*>
*> \param[in] DX
*> \verbatim
*>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>         storage spacing between elements of DX
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup aux_blas
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, linpack, 3/11/78.
*>     modified 3/93 to return if incx .le. 0.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      INTEGER FUNCTION idamax(N,DX,INCX)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER incx,n
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION dx(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION dmax
      INTEGER i,ix
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dabs
*     ..
      idamax = 0
      IF (n.LT.1 .OR. incx.LE.0) RETURN
      idamax = 1
      IF (n.EQ.1) RETURN
      IF (incx.EQ.1) THEN
*
*        code for increment equal to 1
*
         dmax = dabs(dx(1))
         DO i = 2,n
            IF (dabs(dx(i)).GT.dmax) THEN
               idamax = i
               dmax = dabs(dx(i))
            END IF
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         ix = 1
         dmax = dabs(dx(1))
         ix = ix + incx
         DO i = 2,n
            IF (dabs(dx(ix)).GT.dmax) THEN
               idamax = i
               dmax = dabs(dx(ix))
            END IF
            ix = ix + incx
         END DO
      END IF
      RETURN
*
*     End of IDAMAX
*
      END




*> \brief \b ILAENV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,
*                        N4 )
*
*       .. Scalar Arguments ..
*       CHARACTER*( * )    NAME, OPTS
*       INTEGER            ISPEC, N1, N2, N3, N4
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ILAENV returns problem-dependent parameters for the local
*> environment.  See ISPEC for a description of the parameters.
*>
*> In this version, the problem-dependent parameters are contained in
*> the integer array IPARMS in the common block CLAENV and the value
*> with index ISPEC is copied to ILAENV.  This version of ILAENV is
*> to be used in conjunction with XLAENV in TESTING and TIMING.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ISPEC
*> \verbatim
*>          ISPEC is INTEGER
*>          Specifies the parameter to be returned as the value of
*>          ILAENV.
*>          = 1: the optimal blocksize; if this value is 1, an unblocked
*>               algorithm will give the best performance.
*>          = 2: the minimum block size for which the block routine
*>               should be used; if the usable block size is less than
*>               this value, an unblocked routine should be used.
*>          = 3: the crossover point (in a block routine, for N less
*>               than this value, an unblocked routine should be used)
*>          = 4: the number of shifts, used in the nonsymmetric
*>               eigenvalue routines
*>          = 5: the minimum column dimension for blocking to be used;
*>               rectangular blocks must have dimension at least k by m,
*>               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*>          = 6: the crossover point for the SVD (when reducing an m by n
*>               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*>               this value, a QR factorization is used first to reduce
*>               the matrix to a triangular form.)
*>          = 7: the number of processors
*>          = 8: the crossover point for the multishift QR and QZ methods
*>               for nonsymmetric eigenvalue problems.
*>          = 9: maximum size of the subproblems at the bottom of the
*>               computation tree in the divide-and-conquer algorithm
*>          =10: ieee NaN arithmetic can be trusted not to trap
*>          =11: infinity arithmetic can be trusted not to trap
*>
*>          Other specifications (up to 100) can be added later.
*> \endverbatim
*>
*> \param[in] NAME
*> \verbatim
*>          NAME is CHARACTER*(*)
*>          The name of the calling subroutine.
*> \endverbatim
*>
*> \param[in] OPTS
*> \verbatim
*>          OPTS is CHARACTER*(*)
*>          The character options to the subroutine NAME, concatenated
*>          into a single character string.  For example, UPLO = 'U',
*>          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*>          be specified as OPTS = 'UTN'.
*> \endverbatim
*>
*> \param[in] N1
*> \verbatim
*>          N1 is INTEGER
*> \endverbatim
*>
*> \param[in] N2
*> \verbatim
*>          N2 is INTEGER
*> \endverbatim
*>
*> \param[in] N3
*> \verbatim
*>          N3 is INTEGER
*> \endverbatim
*>
*> \param[in] N4
*> \verbatim
*>          N4 is INTEGER
*>
*>          Problem dimensions for the subroutine NAME; these may not all
*>          be required.
*> \endverbatim
*>
*> \return ILAENV
*> \verbatim
*>          ILAENV is INTEGER
*>          >= 0: the value of the parameter specified by ISPEC
*>          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup aux_lin
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The following conventions have been used when calling ILAENV from the
*>  LAPACK routines:
*>  1)  OPTS is a concatenation of all of the character options to
*>      subroutine NAME, in the same order that they appear in the
*>      argument list for NAME, even if they are not used in determining
*>      the value of the parameter specified by ISPEC.
*>  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*>      that they appear in the argument list for NAME.  N1 is used
*>      first, N2 second, and so on, and unused problem dimensions are
*>      passed a value of -1.
*>  3)  The parameter value returned by ILAENV is checked for validity in
*>      the calling subroutine.  For example, ILAENV is used to retrieve
*>      the optimal blocksize for STRTRI as follows:
*>
*>      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*>      IF( NB.LE.1 ) NB = MAX( 1, N )
*> \endverbatim
*>
*  =====================================================================
      INTEGER          FUNCTION ilaenv( ISPEC, NAME, OPTS, N1, N2, N3,
     $                 N4 )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    name, opts
      INTEGER            ispec, n1, n2, n3, n4
*     ..
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          int, min, real
*     ..
*     .. External Functions ..
      INTEGER            ieeeck
      EXTERNAL           ieeeck
*     ..
*     .. Arrays in Common ..
      INTEGER            iparms( 100 )
*     ..
*     .. Common blocks ..
      COMMON             / claenv / iparms
*     ..
*     .. Save statement ..
      SAVE               / claenv /
*     ..
*     .. Executable Statements ..
*
      IF( ispec.GE.1 .AND. ispec.LE.5 ) THEN
*
*        Return a value from the common block.
*
         IF ( name(2:6).EQ.'GEQR ' ) THEN
            IF (n3.EQ.2) THEN
               ilaenv = iparms( 2 )
            ELSE
               ilaenv = iparms( 1 )
            END IF
         ELSE IF ( name(2:6).EQ.'GELQ ' ) THEN
            IF (n3.EQ.2) THEN
               ilaenv = iparms( 2 )
            ELSE
               ilaenv = iparms( 1 )
            END IF
         ELSE
            ilaenv = iparms( ispec )
         END IF
*
      ELSE IF( ispec.EQ.6 ) THEN
*
*        Compute SVD crossover point.
*
         ilaenv = int( real( min( n1, n2 ) )*1.6e0 )
*
      ELSE IF( ispec.GE.7 .AND. ispec.LE.9 ) THEN
*
*        Return a value from the common block.
*
         ilaenv = iparms( ispec )
*
      ELSE IF( ispec.EQ.10 ) THEN
*
*        IEEE NaN arithmetic can be trusted not to trap
*
C        ILAENV = 0
         ilaenv = 1
         IF( ilaenv.EQ.1 ) THEN
            ilaenv = ieeeck( 1, 0.0, 1.0 )
         END IF
*
      ELSE IF( ispec.EQ.11 ) THEN
*
*        Infinity arithmetic can be trusted not to trap
*
C        ILAENV = 0
         ilaenv = 1
         IF( ilaenv.EQ.1 ) THEN
            ilaenv = ieeeck( 0, 0.0, 1.0 )
         END IF
*
      ELSE
*
*        Invalid value for ISPEC
*
         ilaenv = -1
      END IF
*
      RETURN
*
*     End of ILAENV
*
      END
!       INTEGER FUNCTION ilaenv2stage( ISPEC, NAME, OPTS, N1, N2,
!      $                               N3, N4 )
! *     .. Scalar Arguments ..
!       CHARACTER*( * )    name, opts
!       INTEGER            ispec, n1, n2, n3, n4
! *     ..
! *
! *  =====================================================================
! *
! *     .. Local variables ..
!       INTEGER            iispec
! *     .. External Functions ..
!       INTEGER            iparam2stage
!       EXTERNAL           iparam2stage
! *     ..
! *     .. Arrays in Common ..
!       INTEGER            iparms( 100 )
! *     ..
! *     .. Common blocks ..
!       COMMON             / claenv / iparms
! *     ..
! *     .. Save statement ..
!       SAVE               / claenv /
! *     ..
! *     .. Executable Statements ..
! *
!       IF(( ispec.GE.1 ) .AND. (ispec.LE.5)) THEN
! *
! *     1 <= ISPEC <= 5: 2stage eigenvalues SVD routines. 
! *
!          IF( ispec.EQ.1 ) THEN
!              ilaenv2stage = iparms( 1 )
!          ELSE
!              iispec = 16 + ispec
!              ilaenv2stage = iparam2stage( iispec, name, opts,
!      $                                    n1, n2, n3, n4 ) 
!          ENDIF
! *
!       ELSE
! *
! *        Invalid value for ISPEC
! *
!          ilaenv2stage = -1
!       END IF
! *
!       RETURN
!       END




! *> \brief \b IPARAM2STAGE
! *
! *  =========== DOCUMENTATION ===========
! *
! * Online html documentation available at 
! *            http://www.netlib.org/lapack/explore-html/ 
! *
! *> \htmlonly
! *> Download IPARAM2STAGE + dependencies 
! *> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iparam2stage.F"> 
! *> [TGZ]</a> 
! *> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iparam2stage.F"> 
! *> [ZIP]</a>
! *> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iparam2stage.F"> 
! *> [TXT]</a>
! *> \endhtmlonly 
! *
! *  Definition:
! *  ===========
! *
! *       INTEGER FUNCTION IPARAM2STAGE( ISPEC, NAME, OPTS, 
! *                                    NI, NBI, IBI, NXI )
! *       #if defined(_OPENMP)
! *           use omp_lib
! *       #endif
! *       IMPLICIT NONE
! *
! *       .. Scalar Arguments ..
! *       CHARACTER*( * )    NAME, OPTS
! *       INTEGER            ISPEC, NI, NBI, IBI, NXI
! *
! *> \par Purpose:
! *  =============
! *>
! *> \verbatim
! *>
! *>      This program sets problem and machine dependent parameters
! *>      useful for xHETRD_2STAGE, xHETRD_HE2HB, xHETRD_HB2ST,
! *>      xGEBRD_2STAGE, xGEBRD_GE2GB, xGEBRD_GB2BD 
! *>      and related subroutines for eigenvalue problems. 
! *>      It is called whenever ILAENV is called with 17 <= ISPEC <= 21.
! *>      It is called whenever ILAENV2STAGE is called with 1 <= ISPEC <= 5
! *>      with a direct conversion ISPEC + 16.
! *> \endverbatim
! *
! *  Arguments:
! *  ==========
! *
! *> \param[in] ISPEC
! *> \verbatim
! *>          ISPEC is integer scalar
! *>              ISPEC specifies which tunable parameter IPARAM2STAGE should
! *>              return.
! *>
! *>              ISPEC=17: the optimal blocksize nb for the reduction to
! *>                        BAND
! *>
! *>              ISPEC=18: the optimal blocksize ib for the eigenvectors
! *>                        singular vectors update routine
! *>
! *>              ISPEC=19: The length of the array that store the Housholder 
! *>                        representation for the second stage 
! *>                        Band to Tridiagonal or Bidiagonal
! *>
! *>              ISPEC=20: The workspace needed for the routine in input.
! *>
! *>              ISPEC=21: For future release.
! *> \endverbatim
! *>
! *> \param[in] NAME
! *> \verbatim
! *>          NAME is character string
! *>               Name of the calling subroutine
! *> \endverbatim
! *>
! *> \param[in] OPTS
! *> \verbatim
! *>          OPTS is CHARACTER*(*)
! *>          The character options to the subroutine NAME, concatenated
! *>          into a single character string.  For example, UPLO = 'U',
! *>          TRANS = 'T', and DIAG = 'N' for a triangular routine would
! *>          be specified as OPTS = 'UTN'.
! *> \endverbatim
! *>
! *> \param[in] NI
! *> \verbatim
! *>          NI is INTEGER which is the size of the matrix
! *> \endverbatim
! *>
! *> \param[in] NBI
! *> \verbatim
! *>          NBI is INTEGER which is the used in the reduciton, 
! *>          (e.g., the size of the band), needed to compute workspace
! *>          and LHOUS2.
! *> \endverbatim
! *>
! *> \param[in] IBI
! *> \verbatim
! *>          IBI is INTEGER which represent the IB of the reduciton,
! *>          needed to compute workspace and LHOUS2.
! *> \endverbatim
! *>
! *> \param[in] NXI
! *> \verbatim
! *>          NXI is INTEGER needed in the future release.
! *> \endverbatim
! *
! *  Authors:
! *  ========
! *
! *> \author Univ. of Tennessee 
! *> \author Univ. of California Berkeley 
! *> \author Univ. of Colorado Denver 
! *> \author NAG Ltd. 
! *
! *> \ingroup auxOTHERauxiliary
! *
! *> \par Further Details:
! *  =====================
! *>
! *> \verbatim
! *>
! *>  Implemented by Azzam Haidar.
! *>
! *>  All detail are available on technical report, SC11, SC13 papers.
! *>
! *>  Azzam Haidar, Hatem Ltaief, and Jack Dongarra.
! *>  Parallel reduction to condensed forms for symmetric eigenvalue problems
! *>  using aggregated fine-grained and memory-aware kernels. In Proceedings
! *>  of 2011 International Conference for High Performance Computing,
! *>  Networking, Storage and Analysis (SC '11), New York, NY, USA,
! *>  Article 8 , 11 pages.
! *>  http://doi.acm.org/10.1145/2063384.2063394
! *>
! *>  A. Haidar, J. Kurzak, P. Luszczek, 2013.
! *>  An improved parallel singular value algorithm and its implementation 
! *>  for multicore hardware, In Proceedings of 2013 International Conference
! *>  for High Performance Computing, Networking, Storage and Analysis (SC '13).
! *>  Denver, Colorado, USA, 2013.
! *>  Article 90, 12 pages.
! *>  http://doi.acm.org/10.1145/2503210.2503292
! *>
! *>  A. Haidar, R. Solca, S. Tomov, T. Schulthess and J. Dongarra.
! *>  A novel hybrid CPU-GPU generalized eigensolver for electronic structure 
! *>  calculations based on fine-grained memory aware tasks.
! *>  International Journal of High Performance Computing Applications.
! *>  Volume 28 Issue 2, Pages 196-209, May 2014.
! *>  http://hpc.sagepub.com/content/28/2/196 
! *>
! *> \endverbatim
! *>
! *  =====================================================================
!       INTEGER FUNCTION iparam2stage( ISPEC, NAME, OPTS, 
!      $                              NI, NBI, IBI, NXI )
! #if defined(_OPENMP)
!       use omp_lib
! #endif
!       IMPLICIT NONE
! *
! *  -- LAPACK auxiliary routine --
! *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
! *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
! *
! *     .. Scalar Arguments ..
!       CHARACTER*( * )    name, opts
!       INTEGER            ispec, ni, nbi, ibi, nxi
! *
! *  ================================================================
! *     ..
! *     .. Local Scalars ..
!       INTEGER            i, ic, iz, kd, ib, lhous, lwork, nthreads,
!      $                   factoptnb, qroptnb, lqoptnb
!       LOGICAL            rprec, cprec
!       CHARACTER          prec*1, algo*3, stag*5, subnam*12, vect*1
! *     ..
! *     .. Intrinsic Functions ..
!       INTRINSIC          char, ichar, max
! *     ..
! *     .. External Functions ..
!       INTEGER            ilaenv
!       EXTERNAL           ilaenv
! *     ..
! *     .. Executable Statements ..
! *
! *     Invalid value for ISPEC
! *
!       IF( (ispec.LT.17).OR.(ispec.GT.21) ) THEN
!           iparam2stage = -1
!           RETURN
!       ENDIF
! *
! *     Get the number of threads
! *      
!       nthreads = 1
! #if defined(_OPENMP)
! !$OMP PARALLEL 
!       nthreads = omp_get_num_threads()
! !$OMP END PARALLEL
! #endif
! *      WRITE(*,*) 'IPARAM VOICI NTHREADS ISPEC ',NTHREADS, ISPEC
! *
!       IF( ispec .NE. 19 ) THEN
! *
! *        Convert NAME to upper case if the first character is lower case.
! *
!          iparam2stage = -1
!          subnam = name
!          ic = ichar( subnam( 1: 1 ) )
!          iz = ichar( 'Z' )
!          IF( iz.EQ.90 .OR. iz.EQ.122 ) THEN
! *
! *           ASCII character set
! *
!             IF( ic.GE.97 .AND. ic.LE.122 ) THEN
!                subnam( 1: 1 ) = char( ic-32 )
!                DO 100 i = 2, 12
!                   ic = ichar( subnam( i: i ) )
!                   IF( ic.GE.97 .AND. ic.LE.122 )
!      $               subnam( i: i ) = char( ic-32 )
!   100          CONTINUE
!             END IF
! *
!          ELSE IF( iz.EQ.233 .OR. iz.EQ.169 ) THEN
! *
! *           EBCDIC character set
! *
!             IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR.
!      $          ( ic.GE.145 .AND. ic.LE.153 ) .OR.
!      $          ( ic.GE.162 .AND. ic.LE.169 ) ) THEN
!                subnam( 1: 1 ) = char( ic+64 )
!                DO 110 i = 2, 12
!                   ic = ichar( subnam( i: i ) )
!                   IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR.
!      $                ( ic.GE.145 .AND. ic.LE.153 ) .OR.
!      $                ( ic.GE.162 .AND. ic.LE.169 ) )subnam( i:
!      $                i ) = char( ic+64 )
!   110          CONTINUE
!             END IF
! *
!          ELSE IF( iz.EQ.218 .OR. iz.EQ.250 ) THEN
! *
! *           Prime machines:  ASCII+128
! *
!             IF( ic.GE.225 .AND. ic.LE.250 ) THEN
!                subnam( 1: 1 ) = char( ic-32 )
!                DO 120 i = 2, 12
!                  ic = ichar( subnam( i: i ) )
!                  IF( ic.GE.225 .AND. ic.LE.250 )
!      $             subnam( i: i ) = char( ic-32 )
!   120          CONTINUE
!             END IF
!          END IF
! *
!          prec  = subnam( 1: 1 )
!          algo  = subnam( 4: 6 )
!          stag  = subnam( 8:12 )
!          rprec = prec.EQ.'S' .OR. prec.EQ.'D'
!          cprec = prec.EQ.'C' .OR. prec.EQ.'Z'
! *
! *        Invalid value for PRECISION
! *      
!          IF( .NOT.( rprec .OR. cprec ) ) THEN
!              iparam2stage = -1
!              RETURN
!          ENDIF
!       ENDIF
! *      WRITE(*,*),'RPREC,CPREC ',RPREC,CPREC,
! *     $           '   ALGO ',ALGO,'    STAGE ',STAG
! *      
! *
!       IF (( ispec .EQ. 17 ) .OR. ( ispec .EQ. 18 )) THEN 
! *
! *     ISPEC = 17, 18:  block size KD, IB
! *     Could be also dependent from N but for now it
! *     depend only on sequential or parallel
! *
!          IF( nthreads.GT.4 ) THEN
!             IF( cprec ) THEN
!                kd = 128
!                ib = 32
!             ELSE
!                kd = 160
!                ib = 40
!             ENDIF
!          ELSE IF( nthreads.GT.1 ) THEN
!             IF( cprec ) THEN
!                kd = 64
!                ib = 32
!             ELSE
!                kd = 64
!                ib = 32
!             ENDIF
!          ELSE
!             IF( cprec ) THEN
!                kd = 16
!                ib = 16
!             ELSE
!                kd = 32
!                ib = 16
!             ENDIF
!          ENDIF
!          IF( ispec.EQ.17 ) iparam2stage = kd
!          IF( ispec.EQ.18 ) iparam2stage = ib
! *
!       ELSE IF ( ispec .EQ. 19 ) THEN
! *
! *     ISPEC = 19:  
! *     LHOUS length of the Houselholder representation
! *     matrix (V,T) of the second stage. should be >= 1.
! *
! *     Will add the VECT OPTION HERE next release
!          vect  = opts(1:1)
!          IF( vect.EQ.'N' ) THEN
!             lhous = max( 1, 4*ni )
!          ELSE
! *           This is not correct, it need to call the ALGO and the stage2
!             lhous = max( 1, 4*ni ) + ibi
!          ENDIF
!          IF( lhous.GE.0 ) THEN
!             iparam2stage = lhous
!          ELSE
!             iparam2stage = -1
!          ENDIF
! *
!       ELSE IF ( ispec .EQ. 20 ) THEN
! *
! *     ISPEC = 20: (21 for future use)  
! *     LWORK length of the workspace for 
! *     either or both stages for TRD and BRD. should be >= 1.
! *     TRD:
! *     TRD_stage 1: = LT + LW + LS1 + LS2
! *                  = LDT*KD + N*KD + N*MAX(KD,FACTOPTNB) + LDS2*KD 
! *                    where LDT=LDS2=KD
! *                  = N*KD + N*max(KD,FACTOPTNB) + 2*KD*KD
! *     TRD_stage 2: = (2NB+1)*N + KD*NTHREADS
! *     TRD_both   : = max(stage1,stage2) + AB ( AB=(KD+1)*N )
! *                  = N*KD + N*max(KD+1,FACTOPTNB) 
! *                    + max(2*KD*KD, KD*NTHREADS) 
! *                    + (KD+1)*N
!          lwork        = -1
!          subnam(1:1)  = prec
!          subnam(2:6)  = 'GEQRF'
!          qroptnb      = ilaenv( 1, subnam, ' ', ni, nbi, -1, -1 )
!          subnam(2:6)  = 'GELQF'
!          lqoptnb      = ilaenv( 1, subnam, ' ', nbi, ni, -1, -1 )
! *        Could be QR or LQ for TRD and the max for BRD
!          factoptnb    = max(qroptnb, lqoptnb)
!          IF( algo.EQ.'TRD' ) THEN
!             IF( stag.EQ.'2STAG' ) THEN
!                lwork = ni*nbi + ni*max(nbi+1,factoptnb) 
!      $              + max(2*nbi*nbi, nbi*nthreads) 
!      $              + (nbi+1)*ni
!             ELSE IF( (stag.EQ.'HE2HB').OR.(stag.EQ.'SY2SB') ) THEN
!                lwork = ni*nbi + ni*max(nbi,factoptnb) + 2*nbi*nbi
!             ELSE IF( (stag.EQ.'HB2ST').OR.(stag.EQ.'SB2ST') ) THEN
!                lwork = (2*nbi+1)*ni + nbi*nthreads
!             ENDIF
!          ELSE IF( algo.EQ.'BRD' ) THEN
!             IF( stag.EQ.'2STAG' ) THEN
!                lwork = 2*ni*nbi + ni*max(nbi+1,factoptnb) 
!      $              + max(2*nbi*nbi, nbi*nthreads) 
!      $              + (nbi+1)*ni
!             ELSE IF( stag.EQ.'GE2GB' ) THEN
!                lwork = ni*nbi + ni*max(nbi,factoptnb) + 2*nbi*nbi
!             ELSE IF( stag.EQ.'GB2BD' ) THEN
!                lwork = (3*nbi+1)*ni + nbi*nthreads
!             ENDIF
!          ENDIF
!          lwork = max( 1, lwork )
 
!          IF( lwork.GT.0 ) THEN
!             iparam2stage = lwork
!          ELSE
!             iparam2stage = -1
!          ENDIF
! *
!       ELSE IF ( ispec .EQ. 21 ) THEN
! *
! *     ISPEC = 21 for future use 
!          iparam2stage = nxi
!       ENDIF
! *
! *     ==== End of IPARAM2STAGE ====
! *
!       END



*> \brief \b IEEECK
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download IEEECK + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ieeeck.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ieeeck.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ieeeck.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
*
*       .. Scalar Arguments ..
*       INTEGER            ISPEC
*       REAL               ONE, ZERO
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> IEEECK is called from the ILAENV to verify that Infinity and
*> possibly NaN arithmetic is safe (i.e. will not trap).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ISPEC
*> \verbatim
*>          ISPEC is INTEGER
*>          Specifies whether to test just for infinity arithmetic
*>          or whether to test for infinity and NaN arithmetic.
*>          = 0: Verify infinity arithmetic only.
*>          = 1: Verify infinity and NaN arithmetic.
*> \endverbatim
*>
*> \param[in] ZERO
*> \verbatim
*>          ZERO is REAL
*>          Must contain the value 0.0
*>          This is passed to prevent the compiler from optimizing
*>          away this code.
*> \endverbatim
*>
*> \param[in] ONE
*> \verbatim
*>          ONE is REAL
*>          Must contain the value 1.0
*>          This is passed to prevent the compiler from optimizing
*>          away this code.
*>
*>  RETURN VALUE:  INTEGER
*>          = 0:  Arithmetic failed to produce the correct answers
*>          = 1:  Arithmetic produced the correct answers
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup OTHERauxiliary
*
*  =====================================================================
      INTEGER          FUNCTION ieeeck( ISPEC, ZERO, ONE )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            ispec
      REAL               one, zero
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      REAL               nan1, nan2, nan3, nan4, nan5, nan6, neginf,
     $                   negzro, newzro, posinf
*     ..
*     .. Executable Statements ..
      ieeeck = 1
*
      posinf = one / zero
      IF( posinf.LE.one ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      neginf = -one / zero
      IF( neginf.GE.zero ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      negzro = one / ( neginf+one )
      IF( negzro.NE.zero ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      neginf = one / negzro
      IF( neginf.GE.zero ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      newzro = negzro + zero
      IF( newzro.NE.zero ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      posinf = one / newzro
      IF( posinf.LE.one ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      neginf = neginf*posinf
      IF( neginf.GE.zero ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      posinf = posinf*posinf
      IF( posinf.LE.one ) THEN
         ieeeck = 0
         RETURN
      END IF
*
*
*
*
*     Return if we were only asked to check infinity arithmetic
*
      IF( ispec.EQ.0 )
     $   RETURN
*
      nan1 = posinf + neginf
*
      nan2 = posinf / neginf
*
      nan3 = posinf / posinf
*
      nan4 = posinf*zero
*
      nan5 = neginf*negzro
*
      nan6 = nan5*zero
*
      IF( nan1.EQ.nan1 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      IF( nan2.EQ.nan2 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      IF( nan3.EQ.nan3 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      IF( nan4.EQ.nan4 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      IF( nan5.EQ.nan5 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      IF( nan6.EQ.nan6 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      RETURN
      END
