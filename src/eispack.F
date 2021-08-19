c +------------------------------------------------------------------+
c | From Mike Greenfield:                                            |
c | copied on 7/24/93 from the eispack directory                     |
c |   /usr/local/Software/eispack   on mycenae.cchem.berkeley.edu    |
c |                                                                  |
c +------------------------------------------------------------------+
      SUBROUTINE TRED2(NM,N,A,D,E,Z)  
C                                                                       
      INTEGER I,J,K,L,N,II,NM,JP1                                       
      DOUBLE PRECISION A(NM,N),D(N),E(N),Z(NM,N)                        
      DOUBLE PRECISION F,G,H,HH,SCALE                                   
C                                                                       
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,    
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   
C                                                                       
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A              
C     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING               
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.                            
C                                                                       
C     ON INPUT                                                          
C                                                                       
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          
C          DIMENSION STATEMENT.                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX.                                  
C                                                                       
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE          
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.               
C                                                                       
C     ON OUTPUT                                                         
C                                                                       
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.    
C                                                                       
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL         
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.      
C                                                                       
C        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX                
C          PRODUCED IN THE REDUCTION.                                   
C                                                                       
C        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.            
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,    
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY 
C                                                                       
C     THIS VERSION DATED AUGUST 1983.                                   
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      DO 100 I = 1, N                                                   
C                                                                       
         DO 80 J = I, N                                                 
   80    Z(J,I) = A(J,I)                                                
C                                                                       
         D(I) = A(N,I)                                                  
  100 CONTINUE                                                          
C                                                                       
      IF (N .EQ. 1) GO TO 510                                           
C     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........               
      DO 300 II = 2, N                                                  
         I = N + 2 - II                                                 
         L = I - 1                                                      
         H = 0.0D0                                                      
         SCALE = 0.0D0                                                  
         IF (L .LT. 2) GO TO 130                                        
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........       
         DO 120 K = 1, L                                                
  120    SCALE = SCALE + DABS(D(K))                                     
C                                                                       
         IF (SCALE .NE. 0.0D0) GO TO 140                                
  130    E(I) = D(L)                                                    
C                                                                       
         DO 135 J = 1, L                                                
            D(J) = Z(L,J)                                               
            Z(I,J) = 0.0D0                                              
            Z(J,I) = 0.0D0                                              
  135    CONTINUE                                                       
C                                                                       
         GO TO 290                                                      
C                                                                       
  140    DO 150 K = 1, L                                                
            D(K) = D(K) / SCALE                                         
            H = H + D(K) * D(K)                                         
  150    CONTINUE                                                       
C                                                                       
         F = D(L)                                                       
         G = -DSIGN(DSQRT(H),F)                                         
         E(I) = SCALE * G                                               
         H = H - F * G                                                  
         D(L) = F - G                                                   
C     .......... FORM A*U ..........                                    
         DO 170 J = 1, L                                                
  170    E(J) = 0.0D0                                                   
C                                                                       
         DO 240 J = 1, L                                                
            F = D(J)                                                    
            Z(J,I) = F                                                  
            G = E(J) + Z(J,J) * F                                       
            JP1 = J + 1                                                 
            IF (L .LT. JP1) GO TO 220                                   
C                                                                       
            DO 200 K = JP1, L                                           
               G = G + Z(K,J) * D(K)                                    
               E(K) = E(K) + Z(K,J) * F                                 
  200       CONTINUE                                                    
C                                                                       
  220       E(J) = G                                                    
  240    CONTINUE                                                       
C     .......... FORM P ..........                                      
         F = 0.0D0                                                      
C                                                                       
         DO 245 J = 1, L                                                
            E(J) = E(J) / H                                             
            F = F + E(J) * D(J)                                         
  245    CONTINUE                                                       
C                                                                       
         HH = F / (H + H)                                               
C     .......... FORM Q ..........                                      
         DO 250 J = 1, L                                                
  250    E(J) = E(J) - HH * D(J)                                        
C     .......... FORM REDUCED A ..........                              
         DO 280 J = 1, L                                                
            F = D(J)                                                    
            G = E(J)                                                    
C                                                                       
            DO 260 K = J, L                                             
  260       Z(K,J) = Z(K,J) - F * E(K) - G * D(K)                       
C                                                                       
            D(J) = Z(L,J)                                               
            Z(I,J) = 0.0D0                                              
  280    CONTINUE                                                       
C                                                                       
  290    D(I) = H                                                       
  300 CONTINUE                                                          
C     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........     
      DO 500 I = 2, N                                                   
         L = I - 1                                                      
         Z(N,L) = Z(L,L)                                                
         Z(L,L) = 1.0D0                                                 
         H = D(I)                                                       
         IF (H .EQ. 0.0D0) GO TO 380                                    
C                                                                       
         DO 330 K = 1, L                                                
  330    D(K) = Z(K,I) / H                                              
C                                                                       
         DO 360 J = 1, L                                                
            G = 0.0D0                                                   
C                                                                       
            DO 340 K = 1, L                                             
  340       G = G + Z(K,I) * Z(K,J)                                     
C                                                                       
            DO 360 K = 1, L                                             
               Z(K,J) = Z(K,J) - G * D(K)                               
  360    CONTINUE                                                       
C                                                                       
  380    DO 400 K = 1, L                                                
  400    Z(K,I) = 0.0D0                                                 
C                                                                       
  500 CONTINUE                                                          
C                                                                       
  510 DO 520 I = 1, N                                                   
         D(I) = Z(N,I)                                                  
         Z(N,I) = 0.0D0                                                 
  520 CONTINUE                                                          
C                                                                       
      Z(N,N) = 1.0D0                                                    
      E(1) = 0.0D0                                                      
      RETURN                                                            
      END


c +------------------------------------------------------------------+  
c | From Mike Greenfield:                                            |  
c | copied on 7/24/93 from the eispack directory                     |  
c |   /usr/local/Software/eispack   on mycenae.cchem.berkeley.edu    |  
c |                                                                  |  
c +------------------------------------------------------------------+  
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)                                  
C                                                                       
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR                          
      DOUBLE PRECISION D(N),E(N),Z(NM,N)                                
      DOUBLE PRECISION C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2,TST1,TST2,PYTHAG  
C                                                                       
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,     
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND     
C     WILKINSON.                                                        
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).   
C                                                                       
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS            
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.               
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO              
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS                  
C     FULL MATRIX TO TRIDIAGONAL FORM.                                  
C                                                                       
C     ON INPUT                                                          
C                                                                       
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          
C          DIMENSION STATEMENT.                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX.                                  
C                                                                       
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.          
C                                                                       
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX        
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.               
C                                                                       
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE           
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS      
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN        
C          THE IDENTITY MATRIX.                                         
C                                                                       
C      ON OUTPUT                                                        
C                                                                       
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT          
C          UNORDERED FOR INDICES 1,2,...,IERR-1.                        
C                                                                       
C        E HAS BEEN DESTROYED.                                          
C                                                                       
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC           
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,     
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED       
C          EIGENVALUES.  THE EIGENVECTORS ARE STORED AS COLUMN
C          VECTORS.      
C                                                                       
C        IERR IS SET TO                                                 
C          ZERO       FOR NORMAL RETURN,                                
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               
C                     DETERMINED AFTER 30 ITERATIONS.                   
C                                                                       
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .                              
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,    
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY 
C                                                                       
C     THIS VERSION DATED AUGUST 1983.                                   
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      IERR = 0                                                          
      IF (N .EQ. 1) GO TO 1001                                          
C                                                                       
      DO 100 I = 2, N                                                   
  100 E(I-1) = E(I)                                                     
C                                                                       
      F = 0.0D0                                                         
      TST1 = 0.0D0                                                      
      E(N) = 0.0D0                                                      
C                                                                       
      DO 240 L = 1, N                                                   
         J = 0                                                          
         H = DABS(D(L)) + DABS(E(L))                                    
         IF (TST1 .LT. H) TST1 = H                                      
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........         
         DO 110 M = L, N                                                
            TST2 = TST1 + DABS(E(M))                                    
            IF (TST2 .EQ. TST1) GO TO 120                               
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT               
C                THROUGH THE BOTTOM OF THE LOOP ..........              
  110    CONTINUE                                                       
C                                                                       
  120    IF (M .EQ. L) GO TO 220                                        
  130    IF (J .EQ. 30) GO TO 1000                                      
         J = J + 1                                                      
C     .......... FORM SHIFT ..........                                  
         L1 = L + 1                                                     
         L2 = L1 + 1                                                    
         G = D(L)                                                       
         P = (D(L1) - G) / (2.0D0 * E(L))                               
         R = PYTHAG(P,1.0D0)                                            
         D(L) = E(L) / (P + DSIGN(R,P))                                 
         D(L1) = E(L) * (P + DSIGN(R,P))                                
         DL1 = D(L1)                                                    
         H = G - D(L)                                                   
         IF (L2 .GT. N) GO TO 145                                       
C                                                                       
         DO 140 I = L2, N                                               
  140    D(I) = D(I) - H                                                
C                                                                       
  145    F = F + H                                                      
C     .......... QL TRANSFORMATION ..........                           
         P = D(M)                                                       
         C = 1.0D0                                                      
         C2 = C                                                         
         EL1 = E(L1)                                                    
         S = 0.0D0                                                      
         MML = M - L                                                    
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........             
         DO 200 II = 1, MML                                             
            C3 = C2                                                     
            C2 = C                                                      
            S2 = S                                                      
            I = M - II                                                  
            G = C * E(I)                                                
            H = C * P                                                   
            R = PYTHAG(P,E(I))                                          
            E(I+1) = S * R                                              
            S = E(I) / R                                                
            C = P / R                                                   
            P = C * D(I) - S * G                                        
            D(I+1) = H + S * (C * G + S * D(I))                         
C     .......... FORM VECTOR ..........                                 
            DO 180 K = 1, N                                             
               H = Z(K,I+1)                                             
               Z(K,I+1) = S * Z(K,I) + C * H                            
               Z(K,I) = C * Z(K,I) - S * H                              
  180       CONTINUE                                                    
C                                                                       
  200    CONTINUE                                                       
C                                                                       
         P = -S * S2 * C3 * EL1 * E(L) / DL1                            
         E(L) = S * P                                                   
         D(L) = C * P                                                   
         TST2 = TST1 + DABS(E(L))                                       
         IF (TST2 .GT. TST1) GO TO 130                                  
  220    D(L) = D(L) + F                                                
  240 CONTINUE                                                          
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........          
      DO 300 II = 2, N                                                  
         I = II - 1                                                     
         K = I                                                          
         P = D(I)                                                       
C                                                                       
         DO 260 J = II, N                                               
            IF (D(J) .GE. P) GO TO 260                                  
            K = J                                                       
            P = D(J)                                                    
  260    CONTINUE                                                       
C                                                                       
         IF (K .EQ. I) GO TO 300                                        
         D(K) = D(I)                                                    
         D(I) = P                                                       
C                                                                       
         DO 280 J = 1, N                                                
            P = Z(J,I)                                                  
            Z(J,I) = Z(J,K)                                             
            Z(J,K) = P                                                  
  280    CONTINUE                                                       
C                                                                       
  300 CONTINUE                                                          
C                                                                       
      GO TO 1001                                                        
C     .......... SET ERROR -- NO CONVERGENCE TO AN                      
C                EIGENVALUE AFTER 30 ITERATIONS ..........              
 1000 IERR = L                                                          
 1001 RETURN                                                            
      END                                                             

  
c +------------------------------------------------------------------+
c | From Mike Greenfield:                                            |
c | copied on 7/25/93 from the eispack directory                     |
c |   /usr/local/Software/eispack   on mycenae.cchem.berkeley.edu    |
c |                                                                  |
c +------------------------------------------------------------------+

      DOUBLE PRECISION FUNCTION PYTHAG(A,B)                             
      DOUBLE PRECISION A,B                                              
C                                                                       
C     FINDS DSQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW  
C                                                                       
      DOUBLE PRECISION P,R,S,T,U                                        
      P = DMAX1(DABS(A),DABS(B))                                        
      IF (P .EQ. 0.0D0) GO TO 20                                        
      R = (DMIN1(DABS(A),DABS(B))/P)**2                                 
   10 CONTINUE                                                          
         T = 4.0D0 + R                                                  
         IF (T .EQ. 4.0D0) GO TO 20                                     
         S = R/T                                                        
         U = 1.0D0 + 2.0D0*S                                            
         P = U*P                                                        
         R = (S/U)**2 * R                                               
      GO TO 10                                                          
   20 PYTHAG = P                                                        
      RETURN                                                            
      END                                                               






