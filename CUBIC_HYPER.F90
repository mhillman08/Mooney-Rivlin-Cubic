
!THIS STACK OF SUBROUTINES COMPUTES THE SECOND PK STRESS AND ASSOCIATED
!LINEARIZATION TERMS FOR IMPLICIT SOLVERS FOR THE CUBIC MODEL OF MOONEY-RIVLIN
!
!THE CONSTANTS ARE A10, A20, A30, AND BULK, FOR THE MODEL
!
!WRITTEN BY MICHAEL CHARLES HILLMAN CIRCA 2015.


       SUBROUTINE GET_INVAR(F,C,C_INV,JAC,I1,I2,I3,I1_BAR,I2_BAR,I3_BAR)
       !PREPROCESSING: GET INVARIENTS and "BAR" VERSIONS
       
       IMPLICIT NONE
       
       !GLOBAL IN-OUT
       DOUBLE PRECISION::F(3,3),C_INV(3,3),C(3,3),I1,I2,I3,I1_BAR,I2_BAR,I3_BAR
       DOUBLE PRECISION::JAC
       
       !LOCAL
       INTEGER::I,J,K
       
       !Jacobian
       JAC =  F(1,1)*(F(3,3) * F(2,2) - F(3,2) * F(2,3)) - &
            F(2,1)*(F(3,3) * F(1,2) - F(3,2) * F(1,3)) + &
            F(3,1)*(F(2,3) * F(1,2) - F(2,2) * F(1,3))
            
       !right Cauchy-Green deformation tensor 
       !C = F^T * F
       
       C=0.0d0
       DO I=1,3
         DO J = 1,3
           DO K=1,3
           C(I,J) = C(I,J) + F(K,I)*F(K,J)
           END DO
         END DO
       END DO
       
       !inverse of right Cauchy-Green deformation tensor 
       C_INV = C
       CALL inver_3(C_INV,3)
       
       !First invariant
       I1 = C(1,1)+C(2,2)+C(3,3)
       
       !Second invariant
       I2 = C(1,1)*C(2,2)+C(2,2)*C(3,3)+C(1,1)*C(3,3) &
            -C(1,2)*C(2,1)-C(2,3)*C(3,2)-C(1,3)*C(3,1)
       
       !Third invariant
       I3 = C(1,1)*(C(3,3) * C(2,2) - C(3,2) * C(2,3)) - &
            C(2,1)*(C(3,3) * C(1,2) - C(3,2) * C(1,3)) + &
            C(3,1)*(C(2,3) * C(1,2) - C(2,2) * C(1,3))
            
       I1_BAR = I1*I3**(-1.0d0/3.0d0)
       I2_BAR = I2*I3**(-2.0d0/3.0d0)
       I3_BAR = I3
            
       
       RETURN
       
       END SUBROUTINE
	   
	   
       SUBROUTINE CUBIC_HYPER_PRE(A10,A20,A30,I1_BAR,BULK,JAC, & !IN
                                  DELTA,Q1,Q2,T11,T22,T12,P) !OUT
       !PREPROCESSING 
       
       IMPLICIT NONE
       
       !GLOBAL IN-OUT
       DOUBLE PRECISION::A10,A20,A30,I1_BAR,DELTA(3,3),Q1,Q2,T11,T22,T12

       DOUBLE PRECISION::BULK,JAC,P
       
       !LOCAL
       DOUBLE PRECISION::C(3,3),C_INV(3,3),ID(3,3)
       INTEGER::I,J,K
       
       DELTA = 0.0D0
       DO I=1,3
       DELTA(I,I) = 1.0D0
       END DO
       
       Q1 = A10 + 2.0D0*A20*(I1_BAR-3.0D0) + 3.0D0*A30*(I1_BAR-3.0D0)*(I1_BAR-3.0D0)
       Q2 = 0.0D0
       T11 = 2.0D0*A20 + 6.0D0*A30*(I1_BAR-3.0D0)
       T22 = 0.0D0
       T12 = 0.0D0
       
       P = BULK*(JAC-1.0D0)
            
            
       
       RETURN
       
       END SUBROUTINE
       
       






       
       SUBROUTINE CUBIC_HYPER_SIJ(SPK,I1,I2,I3,DELTA,Q1,Q2,C,C_INV,P,JAC)
       !SECOND PK STRESS
       
       IMPLICIT NONE
       
       
       !GLOBAL IN-OUT
       DOUBLE PRECISION:: SPK(3,3)
       INTEGER:: I,J
       DOUBLE PRECISION::I1,I2,I3,DELTA(3,3),Q1,Q2,P,JAC
       DOUBLE PRECISION::C(3,3),C_INV(3,3),ID(3,3)
       
       !LOCAL
       DOUBLE PRECISION:: I3OT,I3TT
       
       I3OT = I3**(-1.0D0/3.0D0)
       I3TT = I3**(-2.0D0/3.0D0)
       
       SPK = 2.0D0*(       Q1*I3OT* (DELTA - 1.0D0/3.0D0*I1*C_INV) &
                          +Q2*I3TT* (I1*DELTA-C-2.0D0/3.0D0*I2*C_INV)) &
                          +P*JAC*C_INV


       END SUBROUTINE
       
       
       
       SUBROUTINE CUBIC_HYPER_DSIJ(DSPK,I1,I1X,I2,I2X,I3,I3X,DELTA,Q1,Q1X,C,CX,C_INV,C_INVX,P,PX,JAC,JACX)
       
       IMPLICIT NONE
       
       
       !GLOBAL IN-OUT
       DOUBLE PRECISION:: DSPK(3,3)
       INTEGER:: I,J
       DOUBLE PRECISION::I1,I1X,I2,I2X,I3,I3X,DELTA(3,3),Q1,Q1X,P,PX,JAC,JACX
       DOUBLE PRECISION::C(3,3),CX(3,3),C_INV(3,3),C_INVX(3,3),ID(3,3)
       
       !LOCAL
       DOUBLE PRECISION:: I3OT,I3TT,I3OTX,I3TTX
       
       I3OT = I3**(-1.0D0/3.0D0)
       I3TT = I3**(-2.0D0/3.0D0)
       
       I3OTX = (-1.0D0/3.0D0)*I3**(-4.0D0/3.0D0)*I3X
       I3TTX = (-2.0D0/3.0D0)*I3**(-5.0D0/3.0D0)*I3X
       
       DSPK = 2.0D0*(   (Q1X*I3OT+Q1*I3OTX)*(DELTA - 1.0D0/3.0D0*I1*C_INV)  +                   &
                         Q1*I3OT*           ( - 1.0D0/3.0D0*I1X*C_INV- 1.0D0/3.0D0*I1*C_INVX)) &
                          +PX*JAC*C_INV  &
                          +P*JACX*C_INV  &
                          +P*JAC*C_INVX 

       RETURN
      
       END SUBROUTINE
       
       
       
               						
        SUBROUTINE CUBIC_HYPER_DIJKL(IHAT,JHAT,KHAT,LHAT, &
									        F,C,C_INV,DELTA, &
									        Q1,Q2,T11,T22,T12, &
									        I1,I2,I3,BULK,JAC,DIJLK)  
       !D_IJKL TERM IN LINEARIZATION

        IMPLICIT NONE

        ! GLOBAL ARRAYS
        INTEGER:: DIM
        INTEGER IHAT,JHAT,KHAT,LHAT
        DOUBLE PRECISION:: F(3,3),C(3,3),C_INV(3,3),DELTA(3,3)
        DOUBLE PRECISION:: Q1,Q2,T11,T22,T12,I1,I2,I3,BULK,JAC
        DOUBLE PRECISION:: DIJLK
        
        ! LOCAL ARRAYS
        INTEGER P, Q
        DOUBLE PRECISION:: CPJQL


        DIJLK = 0.D0

        DO P = 1,3
            DO Q = 1,3
            
                CALL CUBIC_HYPER_CIJKL(P,JHAT,Q,LHAT,C,C_INV,DELTA, &
										        Q1,Q2,T11,T22,T12, &
										        I1,I2,I3,BULK,JAC,CPJQL)
                
                DIJLK = DIJLK + F(IHAT,P)*F(KHAT,Q)*CPJQL

            ENDDO
        ENDDO

        RETURN

        END SUBROUTINE



        SUBROUTINE CUBIC_HYPER_CIJKL(II,JJ,KK,LL,C,C_INV,DELTA, &
									        Q1,Q2,T11,T22,T12, &
									        I1,I2,I3,BULK,JAC,CIJLK) 
       !C_IJKL TERM IN LINEARIZATION
        									
        IMPLICIT NONE

        !GLOBAL ARRAYS
        INTEGER:: DIM
        INTEGER:: II,JJ,KK,LL
        DOUBLE PRECISION:: C(3,3),C_INV(3,3),DELTA(3,3)
        DOUBLE PRECISION:: Q1,Q2,T11,T22,T12,I1,I2,I3,BULK,JAC,CIJLK

        !---------------------------------------								
        ! THIS FUNCTION EVALUATES THE VALUE OF ONE ELEMENT OF THE C TENSOR, I.E 
        ! CIJLK = C(II,JJ,KK,LL)
        ! FOLLOWING THE RKPM FOR LARGE DEFORMATION ANALYSIS PAPER	

        CIJLK = (2.D0/3.D0)*Q1*I3**(-1.D0/3.D0)*( -2.D0*(DELTA(II,JJ)*C_INV(KK,LL) + &
            DELTA(KK,LL)*C_INV(II,JJ)) + & ! C_BAR TERM (DEVIATORIC)
            (1.D0/3.D0)*I1*(2.D0*C_INV(II,JJ)*C_INV(KK,LL) + &
            3.D0*C_INV(II,KK)*C_INV(JJ,LL) + 3.D0*C_INV(II,LL)*C_INV(JJ,KK)) ) + &!1ST LINE
            (4.D0/3.D0)*Q2*I3**(-2.D0/3.D0)*( -2.D0*I1*(DELTA(II,JJ)*C_INV(KK,LL) + &
            DELTA(KK,LL)*C_INV(II,JJ)) + 2.D0*(C_INV(II,KK)*C_INV(JJ,LL) + &
            C_INV(II,LL)*C_INV(JJ,KK)) + & !2ND LINE
            I2* ((4.D0/3.D0)*C_INV(II,JJ)*C_INV(KK,LL) +C_INV(II,KK)*C_INV(JJ,LL) + &
            C_INV(II,LL)*C_INV(KK,JJ)) + (3.D0/2.D0)*(2.D0*DELTA(II,JJ)*DELTA(KK,LL) - &
            DELTA(II,KK)*DELTA(JJ,LL) - DELTA(II,LL)*DELTA(KK,JJ)) )+ & !3RD LINE
            4.D0*T11*I3**(-2.D0/3.D0)*(DELTA(II,JJ)-(1.D0/3.D0)*I1*C_INV(II,JJ))* &
            (DELTA(KK,LL)-(1.D0/3.D0)*I1*C_INV(KK,LL)) +  &
            4.D0*T12*I3**(-1.D0)*(DELTA(II,JJ)-(1.D0/3.D0)*I1*C_INV(II,JJ))* &
            (I1*DELTA(KK,LL) - C(KK,LL) - (2.D0/3.D0)*I2*C_INV(KK,LL)) +  & ! 4TH LINE
            4.D0*T22*I3**(-4.D0/3.D0)*(I1*DELTA(II,JJ) - C(II,JJ) - (2.D0/3.D0)*I2*C_INV(II,JJ))* &
            (I1*DELTA(KK,LL) - C(KK,LL) - (2.D0/3.D0)*I2*C_INV(KK,LL)) + & ! 5TH LINE
            BULK*JAC*(JAC-1.0D0)*(C_INV(II,JJ)*C_INV(KK,LL) - & !C_TILDA TERM (VOLUMETRIC)
            C_INV(II,KK)*C_INV(JJ,LL) - C_INV(II,LL)*C_INV(JJ,KK)) +  &
            BULK * JAC**2.D0*C_INV(II,JJ)*C_INV(KK,LL) 
        	
        RETURN						
        									
        END SUBROUTINE

       
       
       		
        SUBROUTINE CUBIC_HYPER_TIJKL(IHAT,JHAT,KHAT,LHAT,SPK,DELTA,TIJLK)  
        !T_IJKL TERM IN LINEARIZATION

        IMPLICIT NONE

        ! GLOBAL ARRAYS
        INTEGER IHAT,JHAT,KHAT,LHAT
        DOUBLE PRECISION:: SPK(3,3),DELTA(3,3),TIJLK
        
        
        TIJLK = SPK(JHAT,LHAT)*DELTA(IHAT,KHAT)

        RETURN

        END SUBROUTINE
        
        
        
        SUBROUTINE CUBIC_HYPER_TT(DELTA,SPK,TT)
        !MATRIX FORM OF T_IJKL TERM IN LINEARIZATION
        
        ! GLOBAL ARRAYS
        DOUBLE PRECISION:: DELTA(3,3),SPK(3,3),TT(4,4)
        
        ! STRAIN ORDERING:
        ! 11
        ! 22
        ! 12
        ! 21
        
        
        ! ROW 1
        CALL CUBIC_HYPER_TIJKL(1,1,1,1,SPK,DELTA,TT(1,1))  
        CALL CUBIC_HYPER_TIJKL(1,1,2,2,SPK,DELTA,TT(1,2))  
        CALL CUBIC_HYPER_TIJKL(1,1,1,2,SPK,DELTA,TT(1,3))  
        CALL CUBIC_HYPER_TIJKL(1,1,2,1,SPK,DELTA,TT(1,4))  
        
        ! ROW 2
        CALL CUBIC_HYPER_TIJKL(2,2,1,1,SPK,DELTA,TT(2,1))  
        CALL CUBIC_HYPER_TIJKL(2,2,2,2,SPK,DELTA,TT(2,2))  
        CALL CUBIC_HYPER_TIJKL(2,2,1,2,SPK,DELTA,TT(2,3))  
        CALL CUBIC_HYPER_TIJKL(2,2,2,1,SPK,DELTA,TT(2,4))  
        
        ! ROW 3
        CALL CUBIC_HYPER_TIJKL(1,2,1,1,SPK,DELTA,TT(3,1))  
        CALL CUBIC_HYPER_TIJKL(1,2,2,2,SPK,DELTA,TT(3,2))  
        CALL CUBIC_HYPER_TIJKL(1,2,1,2,SPK,DELTA,TT(3,3))  
        CALL CUBIC_HYPER_TIJKL(1,2,2,1,SPK,DELTA,TT(3,4))  
        
        ! ROW 4
        CALL CUBIC_HYPER_TIJKL(2,1,1,1,SPK,DELTA,TT(4,1))  
        CALL CUBIC_HYPER_TIJKL(2,1,2,2,SPK,DELTA,TT(4,2))  
        CALL CUBIC_HYPER_TIJKL(2,1,1,2,SPK,DELTA,TT(4,3))  
        CALL CUBIC_HYPER_TIJKL(2,1,2,1,SPK,DELTA,TT(4,4))  
        
        
        RETURN
        
        END SUBROUTINE
        
        
        
        
        
        
        SUBROUTINE CUBIC_HYPER_TT1(SPK,TT1)
        
        ! GLOBAL ARRAYS
        DOUBLE PRECISION:: SPK(2,2),TT1(4,4)
        
        TT1(1,1) = SPK(1,1)
        TT1(1,2) = 0.0d0
        TT1(1,3) = SPK(1,2)
        TT1(1,4) = 0.0d0
        
        TT1(2,1) = 0.0d0 
        TT1(2,2) = SPK(2,2)
        TT1(2,3) = 0.0d0 
        TT1(2,4) = SPK(2,1)
        
        
        TT1(3,1) = SPK(2,1)
        TT1(3,2) = 0.0d0
        TT1(3,3) = SPK(2,2)
        TT1(3,4) = 0.0d0
        
        TT1(4,1) = 0.0d0 
        TT1(4,2) = SPK(1,2)
        TT1(4,3) = 0.0d0 
        TT1(4,4) = SPK(1,1)
        
        RETURN
        
        END SUBROUTINE
        
        
        
        SUBROUTINE CUBIC_HYPER_TT2(SPK,TT2)
        
        ! GLOBAL ARRAYS
        DOUBLE PRECISION:: SPK(2,2),TT2(4,4)
        
        TT2(1,1) = SPK(1,1)
        TT2(1,2) = 0.0d0
        TT2(1,3) = SPK(1,2)*0.5d0
        TT2(1,4) = SPK(1,2)*0.5d0
        
        TT2(2,2) = SPK(2,2)
        TT2(2,3) = SPK(1,2)*0.5d0
        TT2(2,4) = -SPK(1,2)*0.5d0
        
        TT2(3,3) = (SPK(1,1)+SPK(2,2))*0.25d0
        TT2(3,4) = (SPK(2,2)-SPK(1,1))*0.25d0
        
        TT2(4,4) = (SPK(1,1)+SPK(2,2))*0.25d0
        
        
        !SYM PART
        TT2(2,1) = 0.0d0
        TT2(3,1) = SPK(2,1)*0.5d0
        TT2(3,2) = SPK(2,1)*0.5d0
        TT2(4,1) = SPK(2,1)*0.5d0
        TT2(4,2) = -SPK(2,1)*0.5d0
        TT2(4,3) = (SPK(2,2)-SPK(1,1))*0.25d0
        
        
        
        
        
        RETURN
        
        END SUBROUTINE
        
        
               						

        SUBROUTINE CUBIC_HYPER_DD(FIP,FKQ,C,C_INV,DELTA, &
							Q1,Q2,T11,T22,T12, &
							I1,I2,I3,BULK,JAC, &
							DD) 
        !MATRIX FORM OF D_IJKL TERM IN LINEARIZATION

        IMPLICIT NONE

        ! GLOBAL ARRAYS
        INTEGER IHAT,JHAT,KHAT,LHAT
        DOUBLE PRECISION:: FIP(3,3),FKQ(3,3),C(3,3),C_INV(3,3),DELTA(3,3)
        DOUBLE PRECISION:: Q1,Q2,T11,T22,T12,I1,I2,I3,BULK,JAC
        DOUBLE PRECISION:: DD(4,4)
        
        ! LOCAL ARRAYS
        INTEGER P, Q,I,J,K,L
        DOUBLE PRECISION:: CPJQL
        
        
        ! STRAIN ORDERING:
        ! 11
        ! 22
        ! 12
        ! 21
        
        DD = 0.0d0
		DO P=1,3
		DO Q=1,3	        
		
		I=1
		J=1
		K=1
		L=1
        CALL CUBIC_HYPER_CIJKL(P,J,L,Q, &
									        C,C_INV,DELTA, &
									        Q1,Q2,T11,T22,T12, &
									        I1,I2,I3,BULK,JAC,CPJQL)
		DD(1,1) = DD(1,1) + FIP(I,P)*CPJQL*FKQ(K,Q) 
        
		I=1
		J=1
		K=2
		L=2
        CALL CUBIC_HYPER_CIJKL(P,J,L,Q, &
									        C,C_INV,DELTA, &
									        Q1,Q2,T11,T22,T12, &
									        I1,I2,I3,BULK,JAC,CPJQL) 
		DD(1,2) = DD(1,2) + FIP(I,P)*CPJQL*FKQ(K,Q) 
				
				
		
		
		I=1
		J=1
		K=1
		L=2					        
        CALL CUBIC_HYPER_CIJKL(P,J,L,Q, &
									        C,C_INV,DELTA, &
									        Q1,Q2,T11,T22,T12, &
									        I1,I2,I3,BULK,JAC,CPJQL)  
		DD(1,3) = DD(1,3) + FIP(I,P)*CPJQL*FKQ(K,Q) 
						
						
						
		
		I=1
		J=1
		K=2
		L=1					        
        CALL CUBIC_HYPER_CIJKL(P,J,L,Q, &
									        C,C_INV,DELTA, &
									        Q1,Q2,T11,T22,T12, &
									        I1,I2,I3,BULK,JAC,CPJQL)  
		DD(1,4) = DD(1,4) + FIP(I,P)*CPJQL*FKQ(K,Q) 
		
							
		I=2
		J=2
		K=2
		L=2					        
        CALL CUBIC_HYPER_CIJKL(P,J,L,Q, &
									        C,C_INV,DELTA, &
									        Q1,Q2,T11,T22,T12, &
									        I1,I2,I3,BULK,JAC,CPJQL)  
		DD(2,2) = DD(2,2) + FIP(I,P)*CPJQL*FKQ(K,Q) 
								
		I=2
		J=2
		K=1
		L=2						        
        CALL CUBIC_HYPER_CIJKL(P,J,L,Q, &
									        C,C_INV,DELTA, &
									        Q1,Q2,T11,T22,T12, &
									        I1,I2,I3,BULK,JAC,CPJQL) 
		DD(2,3) = DD(2,3) + FIP(I,P)*CPJQL*FKQ(K,Q) 
					
								
		I=2
		J=2
		K=2
		L=1						        
        CALL CUBIC_HYPER_CIJKL(P,J,L,Q, &
									        C,C_INV,DELTA, &
									        Q1,Q2,T11,T22,T12, &
									        I1,I2,I3,BULK,JAC,CPJQL) 
		DD(2,4) = DD(2,4) + FIP(I,P)*CPJQL*FKQ(K,Q) 
		
		
							
		I=1
		J=2
		K=1
		L=2	
        CALL CUBIC_HYPER_CIJKL(P,J,L,Q, &
									        C,C_INV,DELTA, &
									        Q1,Q2,T11,T22,T12, &
									        I1,I2,I3,BULK,JAC,CPJQL) 
		DD(3,3) = DD(3,3) +FIP(I,P)*CPJQL*FKQ(K,Q) 
						
		I=1
		J=2
		K=2
		L=1	
        CALL CUBIC_HYPER_CIJKL(P,J,L,Q, &
									        C,C_INV,DELTA, &
									        Q1,Q2,T11,T22,T12, &
									        I1,I2,I3,BULK,JAC,CPJQL) 
		DD(3,4) = DD(3,4) +FIP(I,P)*CPJQL*FKQ(K,Q) 
			
					
		I=2
		J=1
		K=2
		L=1	
        CALL CUBIC_HYPER_CIJKL(P,J,L,Q, &
									        C,C_INV,DELTA, &
									        Q1,Q2,T11,T22,T12, &
									        I1,I2,I3,BULK,JAC,CPJQL) 
		DD(4,4) = DD(4,4) +FIP(I,P)*CPJQL*FKQ(K,Q) 
					
		
		END DO		
		
		END DO
			        		        
		DD(2,1) = DD(1,2)	
			        
		DD(3,1) = DD(1,3)	
		DD(3,2) = DD(2,3)    
					        
		DD(4,1) = DD(1,4)	   
		DD(4,2) = DD(2,4) 
		DD(4,3) = DD(3,4)
		
        RETURN

        END SUBROUTINE
        
        
        
        
        
        
        
		
        SUBROUTINE CUBIC_HYPER_CC(F,C,C_INV,DELTA, &
							Q1,Q2,T11,T22,T12, &
							I1,I2,I3,BULK,JAC, &
							DD) 
        !MATRIX FORM OF C_IJKL TERM IN LINEARIZATION

        IMPLICIT NONE

        ! GLOBAL ARRAYS
        INTEGER IHAT,JHAT,KHAT,LHAT
        DOUBLE PRECISION:: F(3,3),C(3,3),C_INV(3,3),DELTA(3,3)
        DOUBLE PRECISION:: Q1,Q2,T11,T22,T12,I1,I2,I3,BULK,JAC
        DOUBLE PRECISION:: DD(3,3)
        
        ! LOCAL ARRAYS
        INTEGER P, Q
        DOUBLE PRECISION:: CPJQL
        
        
									        
        CALL CUBIC_HYPER_CIJKL(1,1,1,1, &
									        C,C_INV,DELTA, &
									        Q1,Q2,T11,T22,T12, &
									        I1,I2,I3,BULK,JAC,DD(1,1))  
        
        CALL CUBIC_HYPER_CIJKL(1,1,2,2, &
									        C,C_INV,DELTA, &
									        Q1,Q2,T11,T22,T12, &
									        I1,I2,I3,BULK,JAC,DD(1,2))  
									        
        CALL CUBIC_HYPER_CIJKL(1,1,1,2, &
									        C,C_INV,DELTA, &
									        Q1,Q2,T11,T22,T12, &
									        I1,I2,I3,BULK,JAC,DD(1,3))  
									        
									        
		DD(2,1) = DD(1,2)
									        
        CALL CUBIC_HYPER_CIJKL(2,2,2,2, &
									        C,C_INV,DELTA, &
									        Q1,Q2,T11,T22,T12, &
									        I1,I2,I3,BULK,JAC,DD(2,2))  
									        
        CALL CUBIC_HYPER_CIJKL(2,2,1,2, &
									        C,C_INV,DELTA, &
									        Q1,Q2,T11,T22,T12, &
									        I1,I2,I3,BULK,JAC,DD(2,3))  
									        		        
		DD(3,1) = DD(1,3)		        
		DD(3,2) = DD(2,3)
		
        CALL CUBIC_HYPER_CIJKL(1,2,1,2, &
									        C,C_INV,DELTA, &
									        Q1,Q2,T11,T22,T12, &
									        I1,I2,I3,BULK,JAC,DD(3,3)) 
									        
									        

        RETURN

        END SUBROUTINE
        
        
        
        
        
        