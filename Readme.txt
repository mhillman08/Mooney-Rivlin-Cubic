
!THE ROUTINES SHOULD BE CALLED AS FOLLOWS

!AT EACH INTEGRATION POINT... COMPUTE DEFORMATION GRADIENT F

!AFTER COMPUTING F, GET THE INVARIENTS AND SOME OTHER LAGRANGIAN QUANTITIES

        CALL GET_INVAR(F,C,C_INV,JAC,I1,I2,I3,I1_BAR,I2_BAR,I3_BAR)

!COMPUTE SOME MORE QUANITTIES NEEDED FOR THE STRESS

        CALL CUBIC_HYPER_PRE(A10,A20,A30,I1_BAR,BULK,JAC, & !IN
                                  DELTA,Q1,Q2,T11,T22,T12,P) !OUT 
								  
!COMPUTE THE SECOND-PK STRESS					  
        
        CALL CUBIC_HYPER_SIJ(SPK,I1,I2,I3,DELTA,Q1,Q2,C,C_INV,P,JAC)

!FOR LINEARIZATION, GET D_IJKL AND T_IJLK MATRIX TERMS. NOTE THE SYMMETRY...

        CALL CUBIC_HYPER_TT1(SPK,TT1)
        
        CALL CUBIC_HYPER_TT2(SPK,TT2)
        
        
        CALL CUBIC_HYPER_DD(F,C,C_INV,DELTA, &
							Q1,Q2,T11,T22,T12, &
							I1,I2,I3,BULK,JAC, &
							DD)
         
        
        DO J = 1,IP !NUMBER OF NEIGHBORS
          JJ=INODE(J) !NEIGHBORS

             !GEOMETRIC "B" MATRICIES
             BMATJ_G1(1,1) = NX(J)
             BMATJ_G1(1,2) = 0.0d0
             
             BMATJ_G1(2,1) = 0.0d0
             BMATJ_G1(2,2) = NY(J)
             
             BMATJ_G1(3,1) = NY(J)
             BMATJ_G1(3,2) = 0.0d0
             
             BMATJ_G1(4,1) = 0.0d0
             BMATJ_G1(4,2) = NX(J)
             
             
             BMATJ_G2(1,1) = NX(J)
             BMATJ_G2(1,2) = 0.0d0
             
             BMATJ_G2(2,1) = 0.0d0
             BMATJ_G2(2,2) = NY(J)
             
             BMATJ_G2(3,1) = NY(J)
             BMATJ_G2(3,2) = NX(J)
             
             BMATJ_G2(4,1) = NY(J)
             BMATJ_G2(4,2) = -NX(J)
             
             
             !MATERIAL "B" MATRICIES
             BMATJ_M(1,1) = NX(J)
             BMATJ_M(1,2) = 0.0d0
             
             BMATJ_M(2,1) = 0.0d0
             BMATJ_M(2,2) = NY(J)
             
             BMATJ_M(3,1) = NY(J)
             BMATJ_M(3,2) = NX(J)
             
             
                MJJ(1)=(jj-1)*2+1
                MJJ(2)=(jj-1)*2+2
                
                
                CALL MAT_MULT(TT1,BMATJ_G1,TTB1,  4, 4, 2)
                CALL MAT_MULT(TT2,BMATJ_G2,TTB2,  4, 4, 2)
                
                CALL MAT_MULT(DD,BMATJ_M,DDB,  3, 3, 2)

! THIS IS THE FIRST PART OF THE LOOP FOR ASSEMBLING THE STIFFNESS... 
! YOU WILL NEED TO FIGURE OUT THE REST. I'M NOT PUBLISHING THE REST
! OF THIS CODE.
!
! I ALSO WILL NOT HELP YOU FIGURE IT OUT SO PLEASE DO NOT BOTHER
! TO CONTACT ME. THE ROUTINES USE STANDARD NOTATIONS IN NONLINEAR
! SOLID MECHANICS TEXTS (e.g., T_IJLK ETC.). GOOD LUCK!