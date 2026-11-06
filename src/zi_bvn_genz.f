* Adapted from pbivnorm 0.6.0, GPL (>= 2).
* Bivariate routine: Alan Genz and Yihong Ge.
* Normal CDF auxiliary: J. L. Schonfelder.
* Original source: github.com/cran/pbivnorm
* Changes: renamed symbols and added ZIBVN wrapper.
      DOUBLE PRECISION FUNCTION ZIBVU( SH, SK, R )
*
*     A function for computing bivariate normal probabilities;
*       developed using 
*         Drezner, Z. and Wesolowsky, G. O. (1989),
*         On the Computation of the Bivariate Normal Integral,
*         J. Stat. Comput. Simul.. 35 pp. 101-107.
*       with extensive modications for double precisions by    
*         Alan Genz and Yihong Ge
*         Department of Mathematics
*         Washington State University
*         Pullman, WA 99164-3113
*         Email : alangenz@wsu.edu
*
* BVN - calculate the probability that X is larger than SH and Y is
*       larger than SK.
*
* Parameters
*
*   SH  REAL, integration limit
*   SK  REAL, integration limit
*   R   REAL, correlation coefficient
*   LG  INTEGER, number of Gauss Rule Points and Weights
*
      DOUBLE PRECISION BVN, SH, SK, R, ZERO, TWOPI 
      INTEGER I, LG, NG
      PARAMETER ( ZERO = 0, TWOPI = 6.283185307179586D0 ) 
      DOUBLE PRECISION X(10,3), W(10,3), AS, A, B, C, D, RS, XS
      DOUBLE PRECISION ZIPHI, SN, ASR, H, K, BS, HS, HK
      SAVE X, W
*     Gauss Legendre Points and Weights, N =  6
      DATA ( W(I,1), X(I,1), I = 1, 3 ) /
     *  0.1713244923791705D+00,-0.9324695142031522D+00,
     *  0.3607615730481384D+00,-0.6612093864662647D+00,
     *  0.4679139345726904D+00,-0.2386191860831970D+00/
*     Gauss Legendre Points and Weights, N = 12
      DATA ( W(I,2), X(I,2), I = 1, 6 ) /
     *  0.4717533638651177D-01,-0.9815606342467191D+00,
     *  0.1069393259953183D+00,-0.9041172563704750D+00,
     *  0.1600783285433464D+00,-0.7699026741943050D+00,
     *  0.2031674267230659D+00,-0.5873179542866171D+00,
     *  0.2334925365383547D+00,-0.3678314989981802D+00,
     *  0.2491470458134029D+00,-0.1252334085114692D+00/
*     Gauss Legendre Points and Weights, N = 20
      DATA ( W(I,3), X(I,3), I = 1, 10 ) /
     *  0.1761400713915212D-01,-0.9931285991850949D+00,
     *  0.4060142980038694D-01,-0.9639719272779138D+00,
     *  0.6267204833410906D-01,-0.9122344282513259D+00,
     *  0.8327674157670475D-01,-0.8391169718222188D+00,
     *  0.1019301198172404D+00,-0.7463319064601508D+00,
     *  0.1181945319615184D+00,-0.6360536807265150D+00,
     *  0.1316886384491766D+00,-0.5108670019508271D+00,
     *  0.1420961093183821D+00,-0.3737060887154196D+00,
     *  0.1491729864726037D+00,-0.2277858511416451D+00,
     *  0.1527533871307259D+00,-0.7652652113349733D-01/
      IF ( ABS(R) .LT. 0.3 ) THEN
         NG = 1
         LG = 3
      ELSE IF ( ABS(R) .LT. 0.75 ) THEN
         NG = 2
         LG = 6
      ELSE 
         NG = 3
         LG = 10
      ENDIF
      H = SH
      K = SK 
      HK = H*K
      BVN = 0
      IF ( ABS(R) .LT. 0.925 ) THEN
         HS = ( H*H + K*K )/2
         ASR = ASIN(R)
         DO I = 1, LG
            SN = SIN(ASR*( X(I,NG)+1 )/2)
            BVN = BVN + W(I,NG)*EXP( ( SN*HK - HS )/( 1 - SN*SN ) )
            SN = SIN(ASR*(-X(I,NG)+1 )/2)
            BVN = BVN + W(I,NG)*EXP( ( SN*HK - HS )/( 1 - SN*SN ) )
         END DO
         BVN = BVN*ASR/(2*TWOPI) + ZIPHI(-H)*ZIPHI(-K) 
      ELSE
         IF ( R .LT. 0 ) THEN
            K = -K
            HK = -HK
         ENDIF
         IF ( ABS(R) .LT. 1 ) THEN
            AS = ( 1 - R )*( 1 + R )
            A = SQRT(AS)
            BS = ( H - K )**2
            C = ( 4 - HK )/8 
            D = ( 12 - HK )/16
            BVN = A*EXP( -(BS/AS + HK)/2 )
     +             *( 1 - C*(BS - AS)*(1 - D*BS/5)/3 + C*D*AS*AS/5 )
            IF ( HK .GT. -160 ) THEN
               B = SQRT(BS)
               BVN = BVN - EXP(-HK/2)*SQRT(TWOPI)*ZIPHI(-B/A)*B
     +                    *( 1 - C*BS*( 1 - D*BS/5 )/3 ) 
            ENDIF
            A = A/2
            DO I = 1, LG
               XS = ( A*(X(I,NG)+1) )**2
               RS = SQRT( 1 - XS )
               BVN = BVN + A*W(I,NG)*
     +              ( EXP( -BS/(2*XS) - HK/(1+RS) )/RS 
     +              - EXP( -(BS/XS+HK)/2 )*( 1 + C*XS*( 1 + D*XS ) ) )
               XS = AS*(-X(I,NG)+1)**2/4
               RS = SQRT( 1 - XS )
               BVN = BVN + A*W(I,NG)*EXP( -(BS/XS + HK)/2 )
     +                    *( EXP( -HK*(1-RS)/(2*(1+RS)) )/RS 
     +                       - ( 1 + C*XS*( 1 + D*XS ) ) )
            END DO
            BVN = -BVN/TWOPI
         ENDIF
         IF ( R .GT. 0 ) BVN =  BVN + ZIPHI( -MAX( H, K ) )
         IF ( R .LT. 0 ) BVN = -BVN + MAX( ZERO, ZIPHI(-H) - ZIPHI(-K) )     
      ENDIF
      ZIBVU = BVN
      END
*
      DOUBLE PRECISION FUNCTION ZIPHI(Z)
*     
*     Normal distribution probabilities accurate to 1d-15.
*     Reference: J.L. Schonfelder, Math Comp 32(1978), pp 1232-1240. 
*     
      INTEGER I, IM
      DOUBLE PRECISION A(0:43), BM, B, BP, P, RTWO, T, XA, Z
      PARAMETER( RTWO = 1.414213562373095048801688724209D0, IM = 24 )
      SAVE A
      DATA ( A(I), I = 0, 43 )/
     &    6.10143081923200417926465815756D-1,
     &   -4.34841272712577471828182820888D-1,
     &    1.76351193643605501125840298123D-1,
     &   -6.0710795609249414860051215825D-2,
     &    1.7712068995694114486147141191D-2,
     &   -4.321119385567293818599864968D-3, 
     &    8.54216676887098678819832055D-4, 
     &   -1.27155090609162742628893940D-4,
     &    1.1248167243671189468847072D-5, 3.13063885421820972630152D-7,      
     &   -2.70988068537762022009086D-7, 3.0737622701407688440959D-8,
     &    2.515620384817622937314D-9, -1.028929921320319127590D-9,
     &    2.9944052119949939363D-11, 2.6051789687266936290D-11,
     &   -2.634839924171969386D-12, -6.43404509890636443D-13,
     &    1.12457401801663447D-13, 1.7281533389986098D-14, 
     &   -4.264101694942375D-15, -5.45371977880191D-16,
     &    1.58697607761671D-16, 2.0899837844334D-17, 
     &   -5.900526869409D-18, -9.41893387554D-19, 2.14977356470D-19, 
     &    4.6660985008D-20, -7.243011862D-21, -2.387966824D-21, 
     &    1.91177535D-22, 1.20482568D-22, -6.72377D-25, -5.747997D-24,
     &   -4.28493D-25, 2.44856D-25, 4.3793D-26, -8.151D-27, -3.089D-27, 
     &    9.3D-29, 1.74D-28, 1.6D-29, -8.0D-30, -2.0D-30 /       
*     
      XA = ABS(Z)/RTWO
      IF ( XA .GT. 100 ) THEN
         P = 0
      ELSE
         T = ( 8*XA - 30 ) / ( 4*XA + 15 )
         BM = 0
         B  = 0
         DO I = IM, 0, -1 
            BP = B
            B  = BM
            BM = T*B - BP  + A(I)
         END DO
         P = EXP( -XA*XA )*( BM - BP )/4
      END IF
      IF ( Z .GT. 0 ) P = 1 - P
      ZIPHI = P
      END
*
      SUBROUTINE ZIBVN(A1,A2,R,P)
      DOUBLE PRECISION A1,A2,R,P,ZIBVU
      P = ZIBVU(-A1,-A2,R)
      END
