!     NewColumbia Distribution - Version 1.0
!     tag: subroutine ../vents2/lavelle/NewCross/new_hs3crt_composite.f - version 6.4.97

      SUBROUTINE HS3CRT (BDXS,BDXF,BDYS,BDYF,BDZS,BDZF,
     1                   LDIMF,MDIMF,F,PERTRB,W)
C
C     PACKAGE HS3CRT, VERSION 1, AUGUST 1985
C
c     made double precision by j.w.lavelle, july, 1992

      integer i,j,k,l,m,n,lp,mp,np,lperod,mdimf,ldimf,ndimf,iw
      real BDXS(MDIMF,*),BDXF(MDIMF,*),BDYS(LDIMF,*),BDYF(LDIMF,*),
     1          BDZS(LDIMF,*),BDZF(LDIMF,*),F(LDIMF,MDIMF,*),W(*)
      real deltax,deltay,deltaz,elmbda,dlxrcp,dlyrcp,dlzrcp,
     1       twdxsq,twdysq,twdzsq,pertrb,s

      DELTAX=W(1)
      L=W(2)
      LP=W(3)
      DELTAY=W(4)
      M=W(5)
      MP=W(6)
      DELTAZ=W(7)
      N=W(8)
      NP=W(9)
      ELMBDA=W(10)
      LPEROD = 0
      IF (LP .GT. 1) LPEROD = 1
      DLXRCP = 1.0d0/DELTAX
      TWDXSQ = 2.0d0/DELTAX**2
      DLYRCP = 1.0d0/DELTAY
      TWDYSQ = 2.0d0/DELTAY**2
      DLZRCP = 1.0d0/DELTAZ
      TWDZSQ = 2.0d0/DELTAZ**2
C
C     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
C
      GO TO (111,102,102,104,104),LP
  102 DO 103 K=1,N
      DO 103 J=1,M
         F(1,J,K) = F(1,J,K)-BDXS(J,K)*TWDXSQ
  103 CONTINUE
      GO TO 106
  104 DO 105 K=1,N
      DO 105 J=1,M
         F(1,J,K) = F(1,J,K)+BDXS(J,K)*DLXRCP
  105 CONTINUE
  106 GO TO (111,107,109,109,107),LP
  107 DO 108 K=1,N
      DO 108 J=1,M
         F(L,J,K) = F(L,J,K)-BDXF(J,K)*TWDXSQ
  108 CONTINUE
      GO TO 111
  109 DO 110 K=1,N
      DO 110 J=1,M
         F(L,J,K) = F(L,J,K)-BDXF(J,K)*DLXRCP
  110 CONTINUE
  111 CONTINUE
C
C     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
C
      GO TO (121,112,112,114,114),MP
  112 DO 113 K=1,N
      DO 113 I=1,L
         F(I,1,K) = F(I,1,K)-BDYS(I,K)*TWDYSQ
  113 CONTINUE
      GO TO 116
  114 DO 115 K=1,N
      DO 115 I=1,L
         F(I,1,K) = F(I,1,K)+BDYS(I,K)*DLYRCP
  115 CONTINUE
  116 GO TO (121,117,119,119,117),MP
  117 DO 118 K=1,N
      DO 118 I=1,L
         F(I,M,K) = F(I,M,K)-BDYF(I,K)*TWDYSQ
  118 CONTINUE
      GO TO 121
  119 DO 120 K=1,N
      DO 120 I=1,L
         F(I,M,K) = F(I,M,K)-BDYF(I,K)*DLYRCP
  120 CONTINUE
  121 CONTINUE
C
C     ENTER BOUNDARY DATA FOR Z-BOUNDARIES.
C
      GO TO (221,212,212,214,214),NP
  212 DO 213 J=1,M
      DO 213 I=1,L
         F(I,J,1) = F(I,J,1)-BDZS(I,J)*TWDZSQ
  213 CONTINUE
      GO TO 216
  214 DO 215 J=1,M
      DO 215 I=1,L
         F(I,J,1) = F(I,J,1)+BDZS(I,J)*DLZRCP
  215 CONTINUE
  216 GO TO (221,217,219,219,217),NP
  217 DO 218 J=1,M
      DO 218 I=1,L
         F(I,J,N) = F(I,J,N)-BDZF(I,J)*TWDZSQ
  218 CONTINUE
      GO TO 221
  219 DO 220 J=1,M
      DO 220 I=1,L
         F(I,J,N) = F(I,J,N)-BDZF(I,J)*DLZRCP
  220 CONTINUE
  221 CONTINUE
      PERTRB = 0.0d0
      IF (ELMBDA.NE.0.0d0) GO TO 133
  126 GO TO (226,133,133,226,133),LP
  226 GO TO (127,133,133,127,133),MP
  127 GO TO (128,133,133,128,133),NP
C
C     FOR SINGULAR PROBLEMS MUST ADJUST DATA TO INSURE THAT A SOLUTION
C     WILL EXIST.
C
  128 CONTINUE
      S = 0.0d0
      DO 230 K=1,N
      DO 130 J=1,M
         DO 129 I=1,L
            S = S+F(I,J,K)
  129    CONTINUE
  130 CONTINUE
  230 CONTINUE
      PERTRB = S/(L*M*N)
      DO 232 K=1,N
      DO 132 J=1,M
         DO 131 I=1,L
            F(I,J,K) = F(I,J,K)-PERTRB
  131    CONTINUE
  132 CONTINUE
  232 CONTINUE
  133 CONTINUE
C
C     ALLOCATE WORK ARRAY W
C
      IW=3*L+11
C
C     SOLVE THE EQUATION
C
      CALL PSTG3D (LDIMF,MDIMF,F,W(IW))
      RETURN
      END

!     ---------------------------------------------------------------------
 
      SUBROUTINE HS3CRI (XS,XF,L,LBDCND,YS,YF,M,MBDCND,ZS,ZF,N,NBDCND,
     1                   ELMBDA,LDIMF,MDIMF,IERROR,W)
C
C     PACKAGE HS3CRT, VERSION 1, AUGUST 1985
C
c     made double precision by j.w.lavelle, july, 1992

      real W(*)
      real xs,xf,ys,yf,zs,zf,elmbda,deltax,deltay,deltaz,c1,c2,
     1   s,st2
      integer ierror,l,lbdcnd,m,mbdcnd,n,nbdcnd,ldimf,mdimf,
     1   lp,mp,np,lperod,ia,ib,ic,iw,i,ierr1
c
C
C     CHECK FOR INVALID PARAMETERS.
C
      IERROR = 0
      IF (XS.GT.XF) IERROR=1
      IF (L.LT.3) IERROR=2
      IF (LBDCND.LT.0 .OR. LBDCND.GT.4) IERROR=3
      IF (YS.GT.YF) IERROR=4
      IF (M.LT.3) IERROR=5
      IF (MBDCND.LT.0 .OR. MBDCND.GT.4) IERROR=6
      IF (ZS.GT.ZF) IERROR=7
      IF (N.LT.3) IERROR=8
      IF (NBDCND.LT.0 .OR. NBDCND.GT.4) IERROR=9
      IF (LDIMF.LT.L) IERROR=11
      IF (MDIMF.LT.M) IERROR=12
      IF (IERROR.NE.0) RETURN
      IF (ELMBDA.GT.0.) IERROR=10
      DELTAX = (XF-XS)/L
      DELTAY = (YF-YS)/M
      DELTAZ = (ZF-ZS)/N
      C1=1./DELTAY**2
      C2=1./DELTAZ**2
      LP = LBDCND+1
      MP = MBDCND+1
      NP = NBDCND+1
      LPEROD=1
      IF (LBDCND.EQ.0) LPEROD=0
C
C     SAVE PARAMETERS FOR HS3CRT IN WORK ARRAY W.
C
      W(1)=DELTAX
      W(2)=L
      W(3)=LP
      W(4)=DELTAY
      W(5)=M
      W(6)=MP
      W(7)=DELTAZ
      W(8)=N
      W(9)=NP
      W(10)=ELMBDA
C
C     DEFINE THE A,B,C COEFFICIENTS IN W-ARRAY.
C
      IA = 11
      IB = IA+L
      IC = IB+L
      IW = IC+L
      S = 1.0d0/DELTAX**2
      ST2 = 2.0d0*S
      DO 101 I=0,L-1
         W(IA+I) = S
         W(IB+I) = -ST2+ELMBDA
         W(IC+I) = S
  101 CONTINUE
      GO TO (111,102,102,104,104),LP
  102 W(IB)= W(IB)-W(IA)
      GO TO 106
  104 W(IB) = W(IB)+W(IA)
  106 GO TO (111,107,109,109,107),LP
  107 W(IC-1) = W(IC-1)-W(IB-1)
      GO TO 111
  109 W(IC-1) = W(IC-1)+W(IB-1)
  111 CONTINUE
      IF (LBDCND .EQ. 0) GO TO 124
      W(IA) = 0.0d0
      W(IC+L-1) = 0.0d0
  124 CONTINUE
C
C     INITIALIZE SOLVER ROUTINE PSTG3D.
C
      CALL PST3DI(LPEROD,L,MBDCND,M,C1,NBDCND,N,C2,W(IA),W(IB),
     1            W(IC),LDIMF,MDIMF,IERR1,W(IW))
      RETURN
      END

!     __________________________________________________________________________________________

      SUBROUTINE PSTG3D (LDIMF,MDIMF,F,W)
C
C     PACKAGE PSTG3D, VERSION 1, AUGUST 1985
C
      integer i,j,k,l,m,n,lp,mp,np,ldimf,mdimf,ia,ic,icfy,icfz,ifctrd,
     1        iwsy,iwsz,ift,iend,lh,lodd
      real       F(LDIMF,MDIMF,*),W(*)

      L=W(1)
      LP=W(2)
      M=W(3)
      MP=W(4)
      N=W(5)
      NP=W(6)
C
C     ALLOCATION OF WORK ARRAY W
C
      IA=7
      IC=IA+2*L
      ICFY=IC+L
      ICFZ=ICFY+4*M
      IFCTRD=ICFZ+4*N
      IWSY=IFCTRD+L*M*N
      IWSZ=IWSY+M+15
      IFT=IWSZ+N+15
C     IEND=IFT+L*M*N
      GO TO (105,114),LP
C
C     REORDER UNKNOWNS WHEN LPEROD = 0.
C
  105 LH = (L+1)/2
      LODD = 1
      IF (2*LH .EQ. L) LODD = 2
      DO 111 J=1,M
         DO 110 K=1,N
           DO 106 I=1,LH-1
               W(I+IFT) = F(LH-I,J,K)-F(LH+I,J,K)
               W(LH+I+IFT) = F(LH-I,J,K)+F(LH+I,J,K)
  106       CONTINUE
            W(LH+IFT) = 2.0d0*F(LH,J,K)
            GO TO (108,107),LODD
  107       W(L+IFT) = 2.0d0*F(L,J,K)
  108       DO 109 I=1,L
               F(I,J,K) = W(I+IFT)
  109       CONTINUE
  110    CONTINUE
  111 CONTINUE
  114 CONTINUE
      IF (LDIMF.EQ.L .AND. MDIMF.EQ.M) GO TO 300
        CALL P3PACK(F,LDIMF,MDIMF,L,M,N,W(IFT))
      CALL PST3D1 (L,M,MP,N,NP,W(IFT),W(ICFY),W(ICFZ),F,
     1             W(IA),W(IC),W(IFCTRD),W(IWSY),W(IWSZ))
        CALL P3UNPK(F,LDIMF,MDIMF,L,M,N,W(IFT))
       GO TO 400
  300 CALL PST3D1 (L,M,MP,N,NP,F,W(ICFY),W(ICFZ),W(IFT),
     *             W(IA),W(IC),W(IFCTRD),W(IWSY),W(IWSZ))
  400 CONTINUE
      GO TO (115,122),LP
  115 DO 121 J=1,M
         DO 120 K=1,N
            DO 116 I=1,LH-1
               W(LH-I+IFT) = 0.5d0*(F(LH+I,J,K)+F(I,J,K))
               W(LH+I+IFT) = 0.5d0*(F(LH+I,J,K)-F(I,J,K))
  116       CONTINUE
            W(LH+IFT) = 0.5d0*F(LH,J,K)
            GO TO (118,117),LODD
  117       W(L+IFT) = 0.5d0*F(L,J,K)
  118       DO 119 I=1,L
               F(I,J,K) = W(I+IFT)
  119       CONTINUE
  120    CONTINUE
  121 CONTINUE
  122 CONTINUE
      RETURN
      END

! --------------------------------------------------------------------------------------------

      SUBROUTINE PST3DI(LPEROD,L,MPEROD,M,C1,NPEROD,N,C2,A,B,C,
     *                  LDIMF,MDIMF,IERROR,W)
C
C     PACKAGE PSTG3D, VERSION 1, AUGUST 1985
C
       integer i,ia,ib,ic,icfy,icfz,ierror,ifctrd,iwsy,iwsz,l,ldimf,
     &         lp,lperod,m,mdimf,mp,mperod,n,np,nperod       
       real A(L),B(L),C(L),W(*),c1,c2


!     type declarations by j.w.lavelle, december, 1995
C
C     SUBROUTINE THAT INITIALIZES SUBROUTINE PSTG3D.
C
C     CONSULT DOCUMENTATION TO PSTG3D FOR DEFINITION OF PARAMETERS.
C
C     CHECK FOR INVALID INPUT.
C
      IERROR=0
      IF (LPEROD.NE.0 .AND. LPEROD.NE.1) IERROR=1
      IF (L.LT.3) IERROR=2
      IF (MPEROD.LT.0 .AND. MPEROD.GT.4) IERROR=3
      IF (M.LT.3) IERROR=4
      IF (NPEROD.LT.0 .AND. NPEROD.GT.4) IERROR=5
      IF (N.LT.3) IERROR=6
      IF (LDIMF.LT.L) IERROR=7
      IF (MDIMF.LT.M) IERROR=8
      IF (LPEROD.EQ.1 .AND. (A(1).NE.0. .OR. C(L).NE.0.)) IERROR=10
      IF (LPEROD.EQ.1) GO TO 200
      DO 100 I=1,L
      IF (A(I).NE.A(1)) GO TO 110
      IF (B(I).NE.B(1)) GO TO 110
      IF (C(I).NE.A(1)) GO TO 110
  100 CONTINUE
      GO TO 200
  110 IERROR=9
  200 IF (IERROR.NE.0) GO TO 400
C
C     ALLOCATE WORK ARRAY W
C
      IA=7
      IB=IA+L
      IC=IB+L
      ICFY=IC+L
      ICFZ=ICFY+4*M
      IFCTRD=ICFZ+4*N
      IWSY=IFCTRD+L*M*N
      IWSZ=IWSY+M+15
C
C     COPY COEFFICIENT ARRAYS A,B, AND C INTO WORK ARRAY.
C
      DO 300 I=0,L-1
      W(IA+I)=A(I+1)
      W(IB+I)=B(I+1)
  300 W(IC+I)=C(I+1)
      LP=LPEROD+1
      MP=MPEROD+1
      NP=NPEROD+1
      CALL PS3DI1(L,LP,M,MP,C1,N,NP,C2,W(IA),W(IB),W(IC),
     *            W(ICFY),W(ICFZ),W(IFCTRD),W(IWSY),W(IWSZ))
C
C     SAVE PARAMETERS FOR SUBROUTINE PSTG3D IN W.
C
      W(1)=L
      W(2)=LP
      W(3)=M
      W(4)=MP
      W(5)=N
      W(6)=NP
  400 CONTINUE
      RETURN
      END

!     _______________________________________________________________________________

      SUBROUTINE P3PACK(F,LDIMF,MDIMF,L,M,N,G)
C
C     PACKAGE PSTG3D, VERSION 1, AUGUST 1985
C
      integer i,if,j,k,l,ldimf,m,mdimf,n
      real F(LDIMF*MDIMF*N),G(L,M,N)

!     type declarations by j.w.lavelle, december, 1995

C
C      THIS SUBROUTINE PACKS THE SUB-ARRAY F(I,J,K), I=1,...,L,
C      J=1,...,M, AND K=1,...,N, INTO THE ARRAY G.
C
      DO 100 K=1,N
      DO 100 J=1,M
      DO 100 I=1,L
  100 G(I,J,K)=F(LDIMF*(MDIMF*(K-1)+J-1)+I)
      IF=LDIMF*(MDIMF*(N-1)+M-1)+LDIMF
      DO 600 K=N,1,-1
      IF (K.EQ.N) GO TO 350
      DO 200 J=MDIMF,M+1,-1
      IF=IF-LDIMF
      DO 200 I=LDIMF,1,-1
  200 F(IF+I)=F(LDIMF*(MDIMF*(K-1)+J-1)+I)
  350 CONTINUE
      DO 400 J=M,1,-1
      IF=IF-(LDIMF-L)
      DO 400 I=LDIMF,L+1,-1
  400 F(IF+I-L)=F(LDIMF*(MDIMF*(K-1)+J-1)+I)
  600 CONTINUE
      RETURN
      END

!     ----------------------------------------------------------------------------------

      SUBROUTINE PST3D1(L,M,MP,N,NP,F,
     1                  CFY,CFZ,FT,A,C,FCTRD,WSAVEY,WSAVEZ)
C
C     PACKAGE PSTG3D, VERSION 1, AUGUST 1985
C
      integer i,ifwrd,j,k,l,m,mp,n, np
      real A(L),C(L),F(L,M,N),CFY(4*M),CFZ(4*N),FT(L,M,N),
     1          FCTRD(M,L,N),WSAVEY(M+15),WSAVEZ(N+15)

!     type declarations by j.w.lavelle, december, 1995

      IFWRD=1
  100 CONTINUE
C
C     TRANSFORM Y
C
      GO TO (205,220,235,240,255),MP
  205 GO TO (210,215),IFWRD
  210 CALL VSRFTF(F,L,M,N,FT,WSAVEY)
      GO TO 300
  215 CALL VSRFTB(F,L,M,N,FT,WSAVEY)
      GO TO 300
  220 GO TO (225,230),IFWRD
  225 CALL VSSINF(F,L,M,N,FT,CFY(1),CFY(M+1),WSAVEY)
      GO TO 300
  230 CALL VSSINB(F,L,M,N,FT,CFY(1),CFY(M+1),WSAVEY)
      GO TO 300
  235 CALL VSSINQ(F,L,M,N,FT,CFY(1),CFY(M+1),CFY(2*M+1),CFY(3*M+1),
     1            WSAVEY)
      GO TO 300
  240 GO TO (245,250),IFWRD
  245 CALL VSCOSF(F,L,M,N,FT,CFY(1),CFY(M+1),WSAVEY)
      GO TO 300
  250 CALL VSCOSB(F,L,M,N,FT,CFY(1),CFY(M+1),WSAVEY)
      GO TO 300
  255 CALL VSCOSQ(F,L,M,N,FT,CFY(1),CFY(M+1),CFY(2*M+1),CFY(3*M+1),
     1            WSAVEY)
  300 CONTINUE
C
C     TRANSFORM IN Z
C
      GO TO (305,320,335,340,355),NP
  305 GO TO (310,315),IFWRD
  310 CALL VSRFTF(F,L,N,M,FT,WSAVEZ)
      GO TO 400
  315 CALL VSRFTB(F,L,N,M,FT,WSAVEZ)
      GO TO 400
  320 GO TO (325,330),IFWRD
  325 CALL VSSINF(F,L,N,M,FT,CFZ(1),CFZ(N+1),WSAVEZ)
      GO TO 400
  330 CALL VSSINB(F,L,N,M,FT,CFZ(1),CFZ(N+1),WSAVEZ)
      GO TO 400
  335 CALL VSSINQ(F,L,N,M,FT,CFZ(1),CFZ(N+1),CFZ(2*N+1),CFZ(3*N+1),
     1            WSAVEZ)
      GO TO 400
  340 GO TO (345,350),IFWRD
  345 CALL VSCOSF(F,L,N,M,FT,CFZ(1),CFZ(N+1),WSAVEZ)
      GO TO 400
  350 CALL VSCOSB(F,L,N,M,FT,CFZ(1),CFZ(N+1),WSAVEZ)
      GO TO 400
  355 CALL VSCOSQ(F,L,N,M,FT,CFZ(1),CFZ(N+1),CFZ(2*N+1),CFZ(3*N+1),
     1            WSAVEZ)
  400 CONTINUE
      IF (IFWRD.NE.1) GO TO 900
C
C     SOLVE TRIDIAGONAL SYSTEMS (PREVIOUSLY FACTORED IN PST3DI) IN X
C
C        FORWARD SUBSTITUTION
C
      DO 410 K=1,N
      DO 410 J=1,M
  410 F(1,J,K)=F(1,J,K)*FCTRD(J,1,K)
      DO 420 I=2,L
      DO 420 K=1,N
      DO 420 J=1,M
  420 F(I,J,K)=(F(I,J,K)-A(I)*F(I-1,J,K))*FCTRD(J,I,K)
C
C        BACKWARD SUBSTITUTION
C
      DO 430 I=L-1,1,-1
      DO 430 K=1,N
      DO 430 J=1,M
  430 F(I,J,K)=F(I,J,K)-C(I)*FCTRD(J,I,K)*F(I+1,J,K)
      IFWRD=2
      GO TO 100
  900 CONTINUE
      RETURN
      END

!     -----------------------------------------------------------------------

      SUBROUTINE P3UNPK(F,LDIMF,MDIMF,L,M,N,G)
C
C     PACKAGE PSTG3D, VERSION 1, AUGUST 1985
C
      integer  i,if,j,k,l,ldimf,m,mdimf,n
      real F(LDIMF*MDIMF*N),G(L,M,N)

!     type declarations by j.w.lavelle, december, 1995
C
C     THIS SUBROUTINE EXPANDS THE ARRAY G OF DIMENSION L X M X N INTO
C     THE ARRAY F OF DIMENSION LDIMF X MDIMF X N.
C
      IF=L*M*N
      DO 600 K=1,N
      DO 300 J=1,M
      DO 200 I=L+1,LDIMF
  200 F(LDIMF*(MDIMF*(K-1)+J-1)+I)=F(IF+I-L)
  300 IF=IF+(LDIMF-L)
      IF (K.EQ.N) GO TO 590
      DO 500 J=M+1,MDIMF
      DO 400 I=1,LDIMF
  400 F(LDIMF*(MDIMF*(K-1)+J-1)+I)=F(IF+I)
  500 IF=IF+LDIMF
  590 CONTINUE
      DO 100 J=1,M
      DO 100 I=1,L
  100 F(LDIMF*(MDIMF*(K-1)+J-1)+I)=G(I,J,K)
  600 CONTINUE
      RETURN
      END

!     ______________________________________________________________________________

      SUBROUTINE PS3DI1(L,LP,M,MP,C1,N,NP,C2,A,B,C,CFY,CFZ,
     1                  FCTRD,WSAVEY,WSAVEZ)
C
C     PACKAGE PSTG3D, VERSION 1, AUGUST 1985
C
      integer i,j,k,l,lh,lodd,lp,m,mp,n,np
      real A(L),B(L),C(L),CFY(4*M),CFZ(4*N),FCTRD(M,L,N),
     1          WSAVEY(M+15),WSAVEZ(N+15),c1,c2,del,den,pi,pimach

!     type declarations by j.w.lavelle, december, 1995

      IF (LP.NE.1) GO TO 30
      LH=(L+1)/2
      LODD=1
      IF (2*LH.EQ.L) LODD=2
      C(LH-1)=0.
      A(LH)=0.
      C(LH)=2.*C(LH)
      GO TO (10,20), LODD
   10 B(LH-1)=B(LH-1)-A(LH-1)
      B(L)=B(L)+A(L)
      GO TO 30
   20 A(L)=C(LH)
   30 CONTINUE
C
C     COMPUTE TRANSFORM ROOTS
C
      PI=PIMACH(1.0)
      DEL=PI/(2*M)
      GO TO (100,102,104,106,104),MP
  100 CONTINUE
      CFY(1)=0.
      CFY(M)=-4.*C1
      DO 101 J=2,M-1,2
      CFY(J)=-4.*C1*SIN(J*DEL)**2
  101 CFY(J+1)=CFY(J)
      GO TO 111
  102 DO 103 J=1,M
  103 CFY(J)=-4.*C1*SIN(J*DEL)**2
      GO TO 111
  104 DO 105 J=1,M
  105 CFY(J)=-4.*C1*SIN((J-.5)*DEL)**2
      GO TO 111
  106 DO 110 J=1,M
  110 CFY(J)=-4.*C1*SIN((J-1)*DEL)**2
  111 CONTINUE
      DEL=PI/(2*N)
      GO TO (113,115,117,119,117),NP
  113 CONTINUE
      CFZ(1)=0.
      CFZ(N)=-4.*C2
      DO 114 K=2,N-1,2
      CFZ(K)=-4.*C2*SIN(K*DEL)**2
  114 CFZ(K+1)=CFZ(K)
      GO TO 122
  115 DO 116 K=1,N
  116 CFZ(K)=-4.*C2*SIN(K*DEL)**2
      GO TO 122
  117 DO 118 K=1,N
  118 CFZ(K)=-4.*C2*SIN((K-.5)*DEL)**2
      GO TO 122
  119 DO 120 K=1,N
  120 CFZ(K)=-4.*C2*SIN((K-1)*DEL)**2
  122 CONTINUE
C
C     FACTOR M*N TRIDIAGONAL SYSTEMS.  FIRST, DO THE POSSIBLY SINGULAR
C     CASE CORRESPONDING TO J = K = 1.
C
C     NOTE THAT FCTRD IS DIMENSIONED AS F(M,L,N)
C
      FCTRD(1,1,1)=1./(B(1)+CFY(1)+CFZ(1))
      DO 130 I=2,L-1
  130 FCTRD(1,I,1)=1./(B(I)+CFY(1)+CFZ(1)-A(I)*C(I-1)*FCTRD(1,I-1,1))
C
C     IF TRIDIAGONAL SYSTEM (...,A(I),B(I),C(I),...) IS SINGULAR THEN
C     FCTRD(1,L,1) IS 1./0.  IF DENOMINATOR IS EXACTLY ZERO, SET
C     FCTRD(1,L,1)=1.
C
      DEN=B(L)+CFY(1)+CFZ(1)-A(L)*C(L-1)*FCTRD(1,L-1,1)
      IF (DEN .EQ. 0.) DEN =1.
      FCTRD(1,L,1)=1./DEN
C
C     FACTOR CASES J=1, K=2,...,N.
C
       DO 142 K=2,N
  142 FCTRD(1,1,K)=1./(B(1)+CFY(1)+CFZ(K))
      DO 144 I=2,L
      DO 144 K=2,N
  144 FCTRD(1,I,K)=1./(B(I)+CFY(1)+CFZ(K)-A(I)*C(I-1)*FCTRD(1,I-1,K))
C
C     FACTOR CASES K=1, J=2,...,M.
C
      DO 146 J=2,M
  146 FCTRD(J,1,1)=1./(B(1)+CFY(J)+CFZ(1))
      DO 148 I=2,L
      DO 148 J=2,M
  148 FCTRD(J,I,1)=1./(B(I)+CFY(J)+CFZ(1)-A(I)*C(I-1)*FCTRD(J,I-1,1))
C
C     FACTOR REMAINING CASES.
C
      DO 150 K=2,N
      DO 150 J=2,M
  150 FCTRD(J,1,K)=1./(B(1)+CFY(J)+CFZ(K))
      DO 152 I=2,L
      DO 152 K=2,N
      DO 152 J=2,M
  152 FCTRD(J,I,K)=1./(B(I)+CFY(J)+CFZ(K)-A(I)*C(I-1)*FCTRD(J,I-1,K))
C
C     INITIALIZE FFT TRANSFORMS AND PRE-PROCESSING COEFFICIENTS IN Y
C
      GO TO (210,220,230,240,250),MP
  210 CALL VSRFTI(M,WSAVEY)
      GO TO 300
  220 CALL VSSINI(M,CFY(1),CFY(M+1),WSAVEY)
      GO TO 300
  230 CALL VSSNQI(M,CFY(1),CFY(M+1),CFY(2*M+1),CFY(3*M+1),WSAVEY)
      GO TO 300
  240 CALL VSCOSI(M,CFY(1),CFY(M+1),WSAVEY)
      GO TO 300
  250 CALL VSCSQI(M,CFY(1),CFY(M+1),CFY(2*M+1),CFY(3*M+1),WSAVEY)
  300 CONTINUE
C
C     INITIALIZE FFT TRANSFORMS AND PRE-PROCESSING COEFFICIENTS IN Z
C
      GO TO (310,320,330,340,350),NP
  310 CALL VSRFTI(N,WSAVEZ)
      GO TO 400
  320 CALL VSSINI(N,CFZ(1),CFZ(N+1),WSAVEZ)
      GO TO 400
  330 CALL VSSNQI(N,CFZ(1),CFZ(N+1),CFZ(2*N+1),CFZ(3*N+1),WSAVEZ)
      GO TO 400
  340 CALL VSCOSI(N,CFZ(1),CFZ(N+1),WSAVEZ)
      GO TO 400
  350 CALL VSCSQI(N,CFZ(1),CFZ(N+1),CFZ(2*N+1),CFZ(3*N+1),WSAVEZ)
  400 CONTINUE
      RETURN
      END

!     -------------------------------------------------------------------------------

      SUBROUTINE VRADB2 (MP,IDO,L1,CC,CH,MDIMC,WA1)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      integer i,ic,ido,idp2,k,l1,m,mdimc,mp
      real  CC(MDIMC,IDO,2,L1)    ,CH(MDIMC,IDO,L1,2),
     1                WA1(IDO)

!     type declarations by j.w.lavelle, december, 1995

      DO 101 K=1,L1
          DO 1001 M=1,MP
         CH(M,1,K,1) = CC(M,1,1,K)+CC(M,IDO,2,K)
         CH(M,1,K,2) = CC(M,1,1,K)-CC(M,IDO,2,K)
 1001     CONTINUE
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
               DO 1002 M=1,MP
            CH(M,I-1,K,1) = CC(M,I-1,1,K)+CC(M,IC-1,2,K)
            CH(M,I,K,1) = CC(M,I,1,K)-CC(M,IC,2,K)
            CH(M,I-1,K,2) = WA1(I-2)*(CC(M,I-1,1,K)-CC(M,IC-1,2,K))
     1      -WA1(I-1)*(CC(M,I,1,K)+CC(M,IC,2,K))
            CH(M,I,K,2) = WA1(I-2)*(CC(M,I,1,K)+CC(M,IC,2,K))+WA1(I-1)
     1      *(CC(M,I-1,1,K)-CC(M,IC-1,2,K))
 1002          CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
          DO 1003 M=1,MP
         CH(M,IDO,K,1) = CC(M,IDO,1,K)+CC(M,IDO,1,K)
         CH(M,IDO,K,2) = -(CC(M,1,2,K)+CC(M,1,2,K))
 1003     CONTINUE
  106 CONTINUE
  107 RETURN
      END

!     -----------------------------------------------------------------------

      SUBROUTINE VRADB5 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2,WA3,WA4)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      integer i,ic,ido,idp2,k,l1,m,mdimc,mp 
      real  CC(MDIMC,IDO,5,L1)    ,CH(MDIMC,IDO,L1,5),
     1             WA1(IDO)     ,WA2(IDO)     ,WA3(IDO)     ,WA4(IDO)
      real arg,ti11,ti12,tr11,tr12,pimach      

!     type declarations by j.w.lavelle, december, 1995

      ARG=2.*PIMACH(1.0)/5.
      TR11=COS(ARG)
      TI11=SIN(ARG)
      TR12=COS(2.*ARG)
      TI12=SIN(2.*ARG)
      DO 101 K=1,L1
      DO 1001 M=1,MP
         CH(M,1,K,1) = CC(M,1,1,K)+2.*CC(M,IDO,2,K)+2.*CC(M,IDO,4,K)
         CH(M,1,K,2) = (CC(M,1,1,K)+TR11*2.*CC(M,IDO,2,K)
     1   +TR12*2.*CC(M,IDO,4,K))-(TI11*2.*CC(M,1,3,K)
     1   +TI12*2.*CC(M,1,5,K))
         CH(M,1,K,3) = (CC(M,1,1,K)+TR12*2.*CC(M,IDO,2,K)
     1   +TR11*2.*CC(M,IDO,4,K))-(TI12*2.*CC(M,1,3,K)
     1   -TI11*2.*CC(M,1,5,K))
         CH(M,1,K,4) = (CC(M,1,1,K)+TR12*2.*CC(M,IDO,2,K)
     1   +TR11*2.*CC(M,IDO,4,K))+(TI12*2.*CC(M,1,3,K)
     1   -TI11*2.*CC(M,1,5,K))
         CH(M,1,K,5) = (CC(M,1,1,K)+TR11*2.*CC(M,IDO,2,K)
     1   +TR12*2.*CC(M,IDO,4,K))+(TI11*2.*CC(M,1,3,K)
     1   +TI12*2.*CC(M,1,5,K))
 1001          CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
      DO 1002 M=1,MP
            CH(M,I-1,K,1) = CC(M,I-1,1,K)+(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +(CC(M,I-1,5,K)+CC(M,IC-1,4,K))
            CH(M,I,K,1) = CC(M,I,1,K)+(CC(M,I,3,K)-CC(M,IC,2,K))
     1      +(CC(M,I,5,K)-CC(M,IC,4,K))
            CH(M,I-1,K,2) = WA1(I-2)*((CC(M,I-1,1,K)+TR11*
     1      (CC(M,I-1,3,K)+CC(M,IC-1,2,K))+TR12
     1      *(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))-(TI11*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))+TI12*(CC(M,I,5,K)+CC(M,IC,4,K))))
     1      -WA1(I-1)*((CC(M,I,1,K)+TR11*(CC(M,I,3,K)-CC(M,IC,2,K))
     1      +TR12*(CC(M,I,5,K)-CC(M,IC,4,K)))+(TI11*(CC(M,I-1,3,K)
     1      -CC(M,IC-1,2,K))+TI12*(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
            CH(M,I,K,2) = WA1(I-2)*((CC(M,I,1,K)+TR11*(CC(M,I,3,K)
     1      -CC(M,IC,2,K))+TR12*(CC(M,I,5,K)-CC(M,IC,4,K)))
     1      +(TI11*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))+TI12
     1      *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))+WA1(I-1)
     1      *((CC(M,I-1,1,K)+TR11*(CC(M,I-1,3,K)
     1      +CC(M,IC-1,2,K))+TR12*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))
     1      -(TI11*(CC(M,I,3,K)+CC(M,IC,2,K))+TI12
     1      *(CC(M,I,5,K)+CC(M,IC,4,K))))
            CH(M,I-1,K,3) = WA2(I-2)
     1      *((CC(M,I-1,1,K)+TR12*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR11*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))-(TI12*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))-TI11*(CC(M,I,5,K)+CC(M,IC,4,K))))
     1     -WA2(I-1)
     1     *((CC(M,I,1,K)+TR12*(CC(M,I,3,K)-
     1      CC(M,IC,2,K))+TR11*(CC(M,I,5,K)-CC(M,IC,4,K)))
     1      +(TI12*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))-TI11
     1      *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
            CH(M,I,K,3) = WA2(I-2)
     1     *((CC(M,I,1,K)+TR12*(CC(M,I,3,K)-
     1      CC(M,IC,2,K))+TR11*(CC(M,I,5,K)-CC(M,IC,4,K)))
     1      +(TI12*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))-TI11
     1      *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
     1      +WA2(I-1)
     1      *((CC(M,I-1,1,K)+TR12*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR11*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))-(TI12*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))-TI11*(CC(M,I,5,K)+CC(M,IC,4,K))))
            CH(M,I-1,K,4) = WA3(I-2)
     1      *((CC(M,I-1,1,K)+TR12*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR11*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))+(TI12*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))-TI11*(CC(M,I,5,K)+CC(M,IC,4,K))))
     1      -WA3(I-1)
     1     *((CC(M,I,1,K)+TR12*(CC(M,I,3,K)-
     1      CC(M,IC,2,K))+TR11*(CC(M,I,5,K)-CC(M,IC,4,K)))
     1      -(TI12*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))-TI11
     1      *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
            CH(M,I,K,4) = WA3(I-2)
     1     *((CC(M,I,1,K)+TR12*(CC(M,I,3,K)-
     1      CC(M,IC,2,K))+TR11*(CC(M,I,5,K)-CC(M,IC,4,K)))
     1      -(TI12*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))-TI11
     1      *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
     1      +WA3(I-1)
     1      *((CC(M,I-1,1,K)+TR12*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR11*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))+(TI12*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))-TI11*(CC(M,I,5,K)+CC(M,IC,4,K))))
            CH(M,I-1,K,5) = WA4(I-2)
     1      *((CC(M,I-1,1,K)+TR11*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR12*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))+(TI11*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))+TI12*(CC(M,I,5,K)+CC(M,IC,4,K))))
     1      -WA4(I-1)
     1      *((CC(M,I,1,K)+TR11*(CC(M,I,3,K)-CC(M,IC,2,K))
     1      +TR12*(CC(M,I,5,K)-CC(M,IC,4,K)))-(TI11*(CC(M,I-1,3,K)
     1      -CC(M,IC-1,2,K))+TI12*(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
            CH(M,I,K,5) = WA4(I-2)
     1      *((CC(M,I,1,K)+TR11*(CC(M,I,3,K)-CC(M,IC,2,K))
     1      +TR12*(CC(M,I,5,K)-CC(M,IC,4,K)))-(TI11*(CC(M,I-1,3,K)
     1      -CC(M,IC-1,2,K))+TI12*(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
     1      +WA4(I-1)
     1      *((CC(M,I-1,1,K)+TR11*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR12*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))+(TI11*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))+TI12*(CC(M,I,5,K)+CC(M,IC,4,K))))
 1002          CONTINUE
  102    CONTINUE
  103 CONTINUE
      RETURN
      END

!     ---------------------------------------------------------------------------

      SUBROUTINE VRADF3 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      integer i,ic,ido,idp2,k,l1,m,mdimc,mp
      real   CH(MDIMC,IDO,3,L1)  ,CC(MDIMC,IDO,L1,3)     ,
     1                WA1(IDO)     ,WA2(IDO)
      real arg,taui,taur,pimach

!     type declarations by j.w.lavelle, december, 1995

      ARG=2.*PIMACH(1.0)/3.
      TAUR=COS(ARG)
      TAUI=SIN(ARG)
      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = CC(M,1,K,1)+(CC(M,1,K,2)+CC(M,1,K,3))
         CH(M,1,3,K) = TAUI*(CC(M,1,K,3)-CC(M,1,K,2))
         CH(M,IDO,2,K) = CC(M,1,K,1)+TAUR*
     1      (CC(M,1,K,2)+CC(M,1,K,3))
 1001    CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DO 1002 M=1,MP
            CH(M,I-1,1,K) = CC(M,I-1,K,1)+((WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))+(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3)))
            CH(M,I,1,K) = CC(M,I,K,1)+((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3)))
            CH(M,I-1,3,K) = (CC(M,I-1,K,1)+TAUR*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))))+(TAUI*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA2(I-2)*
     1       CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3))))
            CH(M,IC-1,2,K) = (CC(M,I-1,K,1)+TAUR*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))))-(TAUI*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA2(I-2)*
     1       CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3))))
            CH(M,I,3,K) = (CC(M,I,K,1)+TAUR*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))))+(TAUI*((WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2))))
            CH(M,IC,2,K) = (TAUI*((WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2))))-(CC(M,I,K,1)+TAUR*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))))
 1002       CONTINUE
  102    CONTINUE
  103 CONTINUE
      RETURN
      END

!     ------------------------------------------------------------------

      SUBROUTINE VRADFG (MP,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,MDIMC,WA)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      integer i,ic,idij,idl1,ido,idp2,ik,ip,ipp2,ipph,is,j,
     &        j2,jc,k,l,l1,lc,m,mdimc,mp,nbd
      real     CH(MDIMC,IDO,L1,IP)   ,CC(MDIMC,IDO,IP,L1)  ,
     1            C1(MDIMC,IDO,L1,IP)    ,C2(MDIMC,IDL1,IP),
     2                CH2(MDIMC,IDL1,IP)           ,WA(IDO)
      real ai1,ai2,ar1,ar1h,ar2,ar2h,arg,dc2,dcp,ds2,dsp,tpi,pimach

!     type declarations by j.w.lavelle, december, 1995

      TPI=2.*PIMACH(1.0)
      ARG = TPI/REAL(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IPPH = (IP+1)/2
      IPP2 = IP+2
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IF (IDO .EQ. 1) GO TO 119
      DO 101 IK=1,IDL1
         DO 1001 M=1,MP
         CH2(M,IK,1) = C2(M,IK,1)
 1001    CONTINUE
  101 CONTINUE
      DO 103 J=2,IP
         DO 102 K=1,L1
            DO 1002 M=1,MP
            CH(M,1,K,J) = C1(M,1,K,J)
 1002       CONTINUE
  102    CONTINUE
  103 CONTINUE
      IF (NBD .GT. L1) GO TO 107
      IS = -IDO
      DO 106 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 105 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 104 K=1,L1
               DO 1004 M=1,MP
               CH(M,I-1,K,J) = WA(IDIJ-1)*C1(M,I-1,K,J)+WA(IDIJ)
     1           *C1(M,I,K,J)
               CH(M,I,K,J) = WA(IDIJ-1)*C1(M,I,K,J)-WA(IDIJ)
     1           *C1(M,I-1,K,J)
 1004          CONTINUE
  104       CONTINUE
  105    CONTINUE
  106 CONTINUE
      GO TO 111
  107 IS = -IDO
      DO 110 J=2,IP
         IS = IS+IDO
         DO 109 K=1,L1
            IDIJ = IS
            DO 108 I=3,IDO,2
               IDIJ = IDIJ+2
               DO 1008 M=1,MP
               CH(M,I-1,K,J) = WA(IDIJ-1)*C1(M,I-1,K,J)+WA(IDIJ)
     1           *C1(M,I,K,J)
               CH(M,I,K,J) = WA(IDIJ-1)*C1(M,I,K,J)-WA(IDIJ)
     1           *C1(M,I-1,K,J)
 1008          CONTINUE
  108       CONTINUE
  109    CONTINUE
  110 CONTINUE
  111 IF (NBD .LT. L1) GO TO 115
      DO 114 J=2,IPPH
         JC = IPP2-J
         DO 113 K=1,L1
            DO 112 I=3,IDO,2
               DO 1012 M=1,MP
               C1(M,I-1,K,J) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
               C1(M,I-1,K,JC) = CH(M,I,K,J)-CH(M,I,K,JC)
               C1(M,I,K,J) = CH(M,I,K,J)+CH(M,I,K,JC)
               C1(M,I,K,JC) = CH(M,I-1,K,JC)-CH(M,I-1,K,J)
 1012          CONTINUE
  112       CONTINUE
  113    CONTINUE
  114 CONTINUE
      GO TO 121
  115 DO 118 J=2,IPPH
         JC = IPP2-J
         DO 117 I=3,IDO,2
            DO 116 K=1,L1
               DO 1016 M=1,MP
               C1(M,I-1,K,J) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
               C1(M,I-1,K,JC) = CH(M,I,K,J)-CH(M,I,K,JC)
               C1(M,I,K,J) = CH(M,I,K,J)+CH(M,I,K,JC)
               C1(M,I,K,JC) = CH(M,I-1,K,JC)-CH(M,I-1,K,J)
 1016          CONTINUE
  116       CONTINUE
  117    CONTINUE
  118 CONTINUE
      GO TO 121
  119 DO 120 IK=1,IDL1
         DO 1020 M=1,MP
         C2(M,IK,1) = CH2(M,IK,1)
 1020    CONTINUE
  120 CONTINUE
  121 DO 123 J=2,IPPH
         JC = IPP2-J
         DO 122 K=1,L1
            DO 1022 M=1,MP
            C1(M,1,K,J) = CH(M,1,K,J)+CH(M,1,K,JC)
            C1(M,1,K,JC) = CH(M,1,K,JC)-CH(M,1,K,J)
 1022       CONTINUE
  122    CONTINUE
  123 CONTINUE
C
      AR1 = 1.
      AI1 = 0.
      DO 127 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 124 IK=1,IDL1
            DO 1024 M=1,MP
            CH2(M,IK,L) = C2(M,IK,1)+AR1*C2(M,IK,2)
            CH2(M,IK,LC) = AI1*C2(M,IK,IP)
 1024       CONTINUE
  124    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 126 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 125 IK=1,IDL1
               DO 1025 M=1,MP
               CH2(M,IK,L) = CH2(M,IK,L)+AR2*C2(M,IK,J)
               CH2(M,IK,LC) = CH2(M,IK,LC)+AI2*C2(M,IK,JC)
 1025          CONTINUE
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      DO 129 J=2,IPPH
         DO 128 IK=1,IDL1
            DO 1028 M=1,MP
            CH2(M,IK,1) = CH2(M,IK,1)+C2(M,IK,J)
 1028       CONTINUE
  128    CONTINUE
  129 CONTINUE
C
      IF (IDO .LT. L1) GO TO 132
      DO 131 K=1,L1
         DO 130 I=1,IDO
            DO 1030 M=1,MP
            CC(M,I,1,K) = CH(M,I,K,1)
 1030       CONTINUE
  130    CONTINUE
  131 CONTINUE
      GO TO 135
  132 DO 134 I=1,IDO
         DO 133 K=1,L1
            DO 1033 M=1,MP
            CC(M,I,1,K) = CH(M,I,K,1)
 1033       CONTINUE
  133    CONTINUE
  134 CONTINUE
  135 DO 137 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 136 K=1,L1
            DO 1036 M=1,MP
            CC(M,IDO,J2-2,K) = CH(M,1,K,J)
            CC(M,1,J2-1,K) = CH(M,1,K,JC)
 1036       CONTINUE
  136    CONTINUE
  137 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IF (NBD .LT. L1) GO TO 141
      DO 140 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 139 K=1,L1
            DO 138 I=3,IDO,2
               IC = IDP2-I
               DO 1038 M=1,MP
               CC(M,I-1,J2-1,K) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
               CC(M,IC-1,J2-2,K) = CH(M,I-1,K,J)-CH(M,I-1,K,JC)
               CC(M,I,J2-1,K) = CH(M,I,K,J)+CH(M,I,K,JC)
               CC(M,IC,J2-2,K) = CH(M,I,K,JC)-CH(M,I,K,J)
 1038          CONTINUE
  138       CONTINUE
  139    CONTINUE
  140 CONTINUE
      RETURN
  141 DO 144 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 143 I=3,IDO,2
            IC = IDP2-I
            DO 142 K=1,L1
               DO 1042 M=1,MP
               CC(M,I-1,J2-1,K) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
               CC(M,IC-1,J2-2,K) = CH(M,I-1,K,J)-CH(M,I-1,K,JC)
               CC(M,I,J2-1,K) = CH(M,I,K,J)+CH(M,I,K,JC)
               CC(M,IC,J2-2,K) = CH(M,I,K,JC)-CH(M,I,K,J)
 1042          CONTINUE
  142       CONTINUE
  143    CONTINUE
  144 CONTINUE
      RETURN
      END

!     -----------------------------------------------------------------------

      SUBROUTINE VSCOSI(N,C1,C2,WSAVE)
C
C     VSFFTPK, VERSION 1, AUGUST 1985
C
      ENTRY VSSINI(N,C1,C2,WSAVE)
       integer i,n
       real C1(N),C2(N),WSAVE(N+15),c,dx,pi,s,pimach

!     type declarations by j.w.lavelle, december, 1995

      PI=PIMACH(1.0)
      DX=PI/(2*N)
C
C     GENERATE A(I)+-B(I)
C
      DO 100 I=1,N
         C=COS((I-1)*DX)
         S=SIN((I-1)*DX)
         C1(I)=.5*(S+C)
  100    C2(I)=.5*(S-C)
C
C     INITIALIZE VRFFTPK ROUTINE
C
      CALL VRFFTI(N,WSAVE)
      RETURN
      END

!     ------------------------------------------------------------

      SUBROUTINE VSRFTB(F,L,M,N,FT,WSAVE)
C
C     VSFFTPK, VERSION 1, AUGUST 1985
C
      integer i,j,k,k1,k2,l,m,n
      real F(L,N*M),FT(L,N,M),WSAVE(M+15)

!     type declarations by j.w.lavelle, december, 1995

C
C     RE-ORDER INPUT
C
      DO 120 K=1,N
      K1=M*(K-1)+1
      K2=M*K
      DO 100 I=1,L
  100 FT(I,K,1)=F(I,K1)
      IF (2*(M/2) .NE. M) GO TO 120
      DO 110 I=1,L
  110 FT(I,K,M)=F(I,K2)
  120 CONTINUE
      DO 200 J=2,M-1,2
      DO 200 K=1,N
      K1=M*(K-1)+J
      K2=K1+1
      DO 200 I=1,L
      FT(I,K,J)=F(I,K1)
  200 FT(I,K,J+1)=F(I,K2)
C
C     REAL, PERIODIC TRANSFORM
C
      CALL VRFFTB(L*N,M,FT,F,L*N,WSAVE)
      DO 300 J=1,M
      DO 300 K=1,N
      K1=N*(J-1)+K
      DO 300 I=1,L
  300 F(I,K1)=FT(I,K,J)
      RETURN
      END

!     -----------------------------------------------------------------------

      SUBROUTINE VSSINB(F,L,M,N,FT,C1,C2,WORK)
C
C     VSFFTPK, VERSION 1, AUGUST 1985
C
      integer i,j,k,k1,k2,l,m,n
      real F(L,M*N),FT(L,N,M),C1(M),C2(M),WORK(M+15)
      real scale

!     type declarations by j.w.lavelle, december, 1995
C
C     PREPROCESSING
C
      DO 100 J=2,M
      DO 100 K=1,N
         K1=M*(K-1)+J-1
         K2=M*K+1-J
         DO 100 I=1,L
  100       FT(I,K,J)=C1(J)*F(I,K1)-C2(J)*F(I,K2)
      DO 200 K=1,N
         K1=M*K
         DO 200 I=1,L
  200       FT(I,K,1) = .5*F(I,K1)
C
C     REAL,PERIODIC ANALYSIS
C
      CALL VRFFTF(L*N,M,FT,F,L*N,WORK)
C
C     POSTPROCESSING
C
      SCALE=SQRT(2.0d0)
      DO 320 K=1,N
         K1=N*(M-1)+K
         DO 300 I=1,L
  300       F(I,K) = SCALE*FT(I,K,1)
      IF (2*(M/2) .NE. M) GO TO 320
      DO 310 I=1,L
  310       F(I,K1) = -SCALE*FT(I,K,M)
  320 CONTINUE
      DO 400 J=2,M-1,2
      DO 400 K=1,N
         K1=N*(J-1)+K
         K2=K1+N
         DO 400 I=1,L
            F(I,K1) = -SCALE*(FT(I,K,J+1)+FT(I,K,J))
  400       F(I,K2)=-SCALE*(FT(I,K,J+1)-FT(I,K,J))
      RETURN
      END

!     ----------------------------------------------------------------

      SUBROUTINE VSSINQ(F,L,M,N,FT,C1,C2,C3,C4,WORK)
C
C     VSFFTPK, VERSION 1, AUGUST 1985
C
      integer i,j,jby2,k,k1,k2,l,m,n
      real F(L,M*N),FT(L,N,M),C1(M),C2(M),C3(M),C4(M),WORK(M+15)

!     type declarations by j.w.lavelle, december, 1995
C
C     PREPROCESSING
C
      DO 120 K=1,N
         K1=M*(K-1)+1
         K2=M*K
         DO 100 I=1,L
  100       FT(I,K,1)=F(I,K1)
      IF (2*(M/2) .NE. M) GO TO 120
      DO 110 I=1,L
  110       FT(I,K,M)=F(I,K2)
  120 CONTINUE
      DO 200 J=2,M-1,2
      JBY2=J/2
      DO 200 K=1,N
         K1=M*(K-1)+J
         K2=K1+1
         DO 200 I=1,L
            FT(I,K,J)   = F(I,K2)*C1(JBY2)-F(I,K1)*C2(JBY2)
  200       FT(I,K,J+1) = -F(I,K2)*C2(JBY2)-F(I,K1)*C1(JBY2)
C
C     REAL,PERIODIC SYNTHESIS
C
      CALL VRFFTB(L*N,M,FT,F,L*N,WORK)
C
C     POSTPROCESSING
C
      DO 300 J=1,M
      DO 300 K=1,N
         K1=N*(J-1)+K
         DO 300 I=1,L
  300       F(I,K1)=C3(J)*FT(I,K,J)-C4(J)*FT(I,K,M+1-J)
      RETURN
      END

!     -----------------------------------------------------------------------

      SUBROUTINE VRADB3 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      integer i,ic,ido,idp2,k,l1,m,mdimc,mp
      real  CC(MDIMC,IDO,3,L1)    ,CH(MDIMC,IDO,L1,3),
     1                WA1(IDO)   ,WA2(IDO)
      real arg,taui,taur,pimach

!     type declarations by j.w.lavelle, december, 1995

      ARG=2.*PIMACH(1.0)/3.
      TAUR=COS(ARG)
      TAUI=SIN(ARG)
      DO 101 K=1,L1
          DO 1001 M=1,MP
         CH(M,1,K,1) = CC(M,1,1,K)+2.*CC(M,IDO,2,K)
         CH(M,1,K,2) = CC(M,1,1,K)+(2.*TAUR)*CC(M,IDO,2,K)
     1   -(2.*TAUI)*CC(M,1,3,K)
         CH(M,1,K,3) = CC(M,1,1,K)+(2.*TAUR)*CC(M,IDO,2,K)
     1   +2.*TAUI*CC(M,1,3,K)
 1001     CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
               DO 1002 M=1,MP
            CH(M,I-1,K,1) = CC(M,I-1,1,K)+(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
            CH(M,I,K,1) = CC(M,I,1,K)+(CC(M,I,3,K)-CC(M,IC,2,K))
            CH(M,I-1,K,2) = WA1(I-2)*
     1 ((CC(M,I-1,1,K)+TAUR*(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))-
     * (TAUI*(CC(M,I,3,K)+CC(M,IC,2,K))))
     2                   -WA1(I-1)*
     3 ((CC(M,I,1,K)+TAUR*(CC(M,I,3,K)-CC(M,IC,2,K)))+
     * (TAUI*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))))
            CH(M,I,K,2) = WA1(I-2)*
     4 ((CC(M,I,1,K)+TAUR*(CC(M,I,3,K)-CC(M,IC,2,K)))+
     8 (TAUI*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))))
     5                  +WA1(I-1)*
     6 ((CC(M,I-1,1,K)+TAUR*(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))-
     8 (TAUI*(CC(M,I,3,K)+CC(M,IC,2,K))))
              CH(M,I-1,K,3) = WA2(I-2)*
     7 ((CC(M,I-1,1,K)+TAUR*(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))+
     8 (TAUI*(CC(M,I,3,K)+CC(M,IC,2,K))))
     8                      -WA2(I-1)*
     9 ((CC(M,I,1,K)+TAUR*(CC(M,I,3,K)-CC(M,IC,2,K)))-
     8 (TAUI*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))))
            CH(M,I,K,3) = WA2(I-2)*
     1 ((CC(M,I,1,K)+TAUR*(CC(M,I,3,K)-CC(M,IC,2,K)))-
     8 (TAUI*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))))
     2                 +WA2(I-1)*
     3 ((CC(M,I-1,1,K)+TAUR*(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))+
     8 (TAUI*(CC(M,I,3,K)+CC(M,IC,2,K))))
 1002          CONTINUE
  102    CONTINUE
  103 CONTINUE
      RETURN
      END

!     --------------------------------------------------------------------------

      SUBROUTINE VRADBG (MP,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
     *                 MDIMC,WA)
      integer i,ic,idij,idl1,ido,idp2,ik,ip,ipp2,ipph,is,j,
     &        j2,jc,k,l,l1,lc,m,mdimc,mp,nbd
      real    CH(MDIMC,IDO,L1,IP)    ,CC(MDIMC,IDO,IP,L1) ,
     1           C1(MDIMC,IDO,L1,IP)     ,C2(MDIMC,IDL1,IP),
     2                CH2(MDIMC,IDL1,IP)       ,WA(IDO)
      real ai1,ai2,ar1,ar1h,ar2,ar2h,arg,dc2,dcp,ds2,dsp,tpi,pimach

!     type declarations by j.w.lavelle, december, 1995

      TPI=2.*PIMACH(1.0)
      ARG = TPI/REAL(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IF (IDO .LT. L1) GO TO 103
      DO 102 K=1,L1
         DO 101 I=1,IDO
            DO 1001 M=1,MP
            CH(M,I,K,1) = CC(M,I,1,K)
 1001       CONTINUE
  101    CONTINUE
  102 CONTINUE
      GO TO 106
  103 DO 105 I=1,IDO
         DO 104 K=1,L1
            DO 1004 M=1,MP
            CH(M,I,K,1) = CC(M,I,1,K)
 1004       CONTINUE
  104    CONTINUE
  105 CONTINUE
  106 DO 108 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 107 K=1,L1
            DO 1007 M=1,MP
            CH(M,1,K,J) = CC(M,IDO,J2-2,K)+CC(M,IDO,J2-2,K)
            CH(M,1,K,JC) = CC(M,1,J2-1,K)+CC(M,1,J2-1,K)
 1007       CONTINUE
  107    CONTINUE
  108 CONTINUE
      IF (IDO .EQ. 1) GO TO 116
      IF (NBD .LT. L1) GO TO 112
      DO 111 J=2,IPPH
         JC = IPP2-J
         DO 110 K=1,L1
            DO 109 I=3,IDO,2
               IC = IDP2-I
               DO 1009 M=1,MP
               CH(M,I-1,K,J) = CC(M,I-1,2*J-1,K)+CC(M,IC-1,2*J-2,K)
               CH(M,I-1,K,JC) = CC(M,I-1,2*J-1,K)-CC(M,IC-1,2*J-2,K)
               CH(M,I,K,J) = CC(M,I,2*J-1,K)-CC(M,IC,2*J-2,K)
               CH(M,I,K,JC) = CC(M,I,2*J-1,K)+CC(M,IC,2*J-2,K)
 1009          CONTINUE
  109       CONTINUE
  110    CONTINUE
  111 CONTINUE
      GO TO 116
  112 DO 115 J=2,IPPH
         JC = IPP2-J
         DO 114 I=3,IDO,2
            IC = IDP2-I
            DO 113 K=1,L1
               DO 1013 M=1,MP
               CH(M,I-1,K,J) = CC(M,I-1,2*J-1,K)+CC(M,IC-1,2*J-2,K)
               CH(M,I-1,K,JC) = CC(M,I-1,2*J-1,K)-CC(M,IC-1,2*J-2,K)
               CH(M,I,K,J) = CC(M,I,2*J-1,K)-CC(M,IC,2*J-2,K)
               CH(M,I,K,JC) = CC(M,I,2*J-1,K)+CC(M,IC,2*J-2,K)
 1013          CONTINUE
  113       CONTINUE
  114    CONTINUE
  115 CONTINUE
  116 AR1 = 1.
      AI1 = 0.
      DO 120 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 117 IK=1,IDL1
            DO 1017 M=1,MP
            C2(M,IK,L) = CH2(M,IK,1)+AR1*CH2(M,IK,2)
            C2(M,IK,LC) = AI1*CH2(M,IK,IP)
 1017       CONTINUE
  117    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 119 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 118 IK=1,IDL1
               DO 1018 M=1,MP
               C2(M,IK,L) = C2(M,IK,L)+AR2*CH2(M,IK,J)
               C2(M,IK,LC) = C2(M,IK,LC)+AI2*CH2(M,IK,JC)
 1018          CONTINUE
  118       CONTINUE
  119    CONTINUE
  120 CONTINUE
      DO 122 J=2,IPPH
         DO 121 IK=1,IDL1
            DO 1021 M=1,MP
            CH2(M,IK,1) = CH2(M,IK,1)+CH2(M,IK,J)
 1021       CONTINUE
  121    CONTINUE
  122 CONTINUE
      DO 124 J=2,IPPH
         JC = IPP2-J
         DO 123 K=1,L1
            DO 1023 M=1,MP
            CH(M,1,K,J) = C1(M,1,K,J)-C1(M,1,K,JC)
            CH(M,1,K,JC) = C1(M,1,K,J)+C1(M,1,K,JC)
 1023       CONTINUE
  123    CONTINUE
  124 CONTINUE
      IF (IDO .EQ. 1) GO TO 132
      IF (NBD .LT. L1) GO TO 128
      DO 127 J=2,IPPH
         JC = IPP2-J
         DO 126 K=1,L1
            DO 125 I=3,IDO,2
               DO 1025 M=1,MP
               CH(M,I-1,K,J) = C1(M,I-1,K,J)-C1(M,I,K,JC)
               CH(M,I-1,K,JC) = C1(M,I-1,K,J)+C1(M,I,K,JC)
               CH(M,I,K,J) = C1(M,I,K,J)+C1(M,I-1,K,JC)
               CH(M,I,K,JC) = C1(M,I,K,J)-C1(M,I-1,K,JC)
 1025          CONTINUE
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      GO TO 132
  128 DO 131 J=2,IPPH
         JC = IPP2-J
         DO 130 I=3,IDO,2
            DO 129 K=1,L1
               DO 1029 M=1,MP
               CH(M,I-1,K,J) = C1(M,I-1,K,J)-C1(M,I,K,JC)
               CH(M,I-1,K,JC) = C1(M,I-1,K,J)+C1(M,I,K,JC)
               CH(M,I,K,J) = C1(M,I,K,J)+C1(M,I-1,K,JC)
               CH(M,I,K,JC) = C1(M,I,K,J)-C1(M,I-1,K,JC)
 1029          CONTINUE
  129       CONTINUE
  130    CONTINUE
  131 CONTINUE
  132 CONTINUE
      IF (IDO .EQ. 1) RETURN
      DO 133 IK=1,IDL1
         DO 1033 M=1,MP
         C2(M,IK,1) = CH2(M,IK,1)
 1033    CONTINUE
  133 CONTINUE
      DO 135 J=2,IP
         DO 134 K=1,L1
            DO 1034 M=1,MP
            C1(M,1,K,J) = CH(M,1,K,J)
 1034       CONTINUE
  134    CONTINUE
  135 CONTINUE
      IF (NBD .GT. L1) GO TO 139
      IS = -IDO
      DO 138 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 137 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 136 K=1,L1
               DO 1036 M=1,MP
               C1(M,I-1,K,J) = WA(IDIJ-1)*CH(M,I-1,K,J)-WA(IDIJ)*
     1          CH(M,I,K,J)
               C1(M,I,K,J) = WA(IDIJ-1)*CH(M,I,K,J)+WA(IDIJ)*
     1          CH(M,I-1,K,J)
 1036          CONTINUE
  136       CONTINUE
  137    CONTINUE
  138 CONTINUE
      GO TO 143
  139 IS = -IDO
      DO 142 J=2,IP
         IS = IS+IDO
         DO 141 K=1,L1
            IDIJ = IS
            DO 140 I=3,IDO,2
               IDIJ = IDIJ+2
               DO 1040 M=1,MP
               C1(M,I-1,K,J) = WA(IDIJ-1)*CH(M,I-1,K,J)-WA(IDIJ)*
     1          CH(M,I,K,J)
               C1(M,I,K,J) = WA(IDIJ-1)*CH(M,I,K,J)+WA(IDIJ)*
     1          CH(M,I-1,K,J)
 1040          CONTINUE
  140       CONTINUE
  141    CONTINUE
  142 CONTINUE
  143 RETURN
      END

!     ---------------------------------------------------------------

      SUBROUTINE VRADF4 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2,WA3)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      integer i,ic,ido,idp2,k,l1,m,mdimc,mp
      real    CC(MDIMC,IDO,L1,4)   ,CH(MDIMC,IDO,4,L1)     ,
     1                WA1(IDO)     ,WA2(IDO)     ,WA3(IDO)
      real hsqt2

!     type declarations by j.w.lavelle, december, 1995

      HSQT2=SQRT(2.)/2.
      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = (CC(M,1,K,2)+CC(M,1,K,4))
     1      +(CC(M,1,K,1)+CC(M,1,K,3))
         CH(M,IDO,4,K) = (CC(M,1,K,1)+CC(M,1,K,3))
     1      -(CC(M,1,K,2)+CC(M,1,K,4))
         CH(M,IDO,2,K) = CC(M,1,K,1)-CC(M,1,K,3)
         CH(M,1,3,K) = CC(M,1,K,4)-CC(M,1,K,2)
 1001    CONTINUE
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            DO 1003 M=1,MP
            CH(M,I-1,1,K) = ((WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2))+(WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4)))+(CC(M,I-1,K,1)+(WA2(I-2)*CC(M,I-1,K,3)+
     1       WA2(I-1)*CC(M,I,K,3)))
            CH(M,IC-1,4,K) = (CC(M,I-1,K,1)+(WA2(I-2)*CC(M,I-1,K,3)+
     1       WA2(I-1)*CC(M,I,K,3)))-((WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))+(WA3(I-2)*CC(M,I-1,K,4)+
     1       WA3(I-1)*CC(M,I,K,4)))
            CH(M,I,1,K) = ((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))+(CC(M,I,K,1)+(WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3)))
            CH(M,IC,4,K) = ((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))-(CC(M,I,K,1)+(WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3)))
            CH(M,I-1,3,K) = ((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))+(CC(M,I-1,K,1)-(WA2(I-2)*CC(M,I-1,K,3)+
     1       WA2(I-1)*CC(M,I,K,3)))
            CH(M,IC-1,2,K) = (CC(M,I-1,K,1)-(WA2(I-2)*CC(M,I-1,K,3)+
     1       WA2(I-1)*CC(M,I,K,3)))-((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))
            CH(M,I,3,K) = ((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))+(CC(M,I,K,1)-(WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3)))
            CH(M,IC,2,K) = ((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))-(CC(M,I,K,1)-(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3)))
 1003       CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
         DO 1006 M=1,MP
            CH(M,IDO,1,K) = (HSQT2*(CC(M,IDO,K,2)-CC(M,IDO,K,4)))+
     1       CC(M,IDO,K,1)
            CH(M,IDO,3,K) = CC(M,IDO,K,1)-(HSQT2*(CC(M,IDO,K,2)-
     1       CC(M,IDO,K,4)))
            CH(M,1,2,K) = (-HSQT2*(CC(M,IDO,K,2)+CC(M,IDO,K,4)))-
     1       CC(M,IDO,K,3)
            CH(M,1,4,K) = (-HSQT2*(CC(M,IDO,K,2)+CC(M,IDO,K,4)))+
     1       CC(M,IDO,K,3)
 1006    CONTINUE
  106 CONTINUE
  107 RETURN
      END

!     --------------------------------------------------------------------------------------

      SUBROUTINE VRFFTB(M,N,R,RT,MDIMR,WSAVE)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      integer m,mdimr,n
      real     R(MDIMR,N),RT(MDIMR,N),WSAVE(N+15)

!     type declarations by j.w.lavelle, december, 1995

      IF (N .EQ. 1) RETURN
      CALL VRFTB1 (M,N,R,RT,MDIMR,WSAVE(1),WSAVE(N+1))
      RETURN
      END

!     ----------------------------------------------------------------

      SUBROUTINE VRFTB1 (M,N,C,CH,MDIMC,WA,FAC)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      integer i,idl1,ido,ip,iw,ix2,ix3,ix4,j,k1,l1,l2,m,mdimc,n,na,nf
      real       CH(MDIMC,N), C(MDIMC,N), WA(N) ,FAC(15)
      real scale
!     type declarations by j.w.lavelle, december, 1995

      NF = FAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = FAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDL1 = IDO*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL VRADB4 (M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL VRADB4 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL VRADB2 (M,IDO,L1,C,CH,MDIMC,WA(IW))
         GO TO 105
  104    CALL VRADB2 (M,IDO,L1,CH,C,MDIMC,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 107
         CALL VRADB3 (M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2))
         GO TO 108
  107    CALL VRADB3 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 110
      CALL VRADB5 (M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110 CALL VRADB5 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL VRADBG (M,IDO,IP,L1,IDL1,C,C,C,CH,CH,MDIMC,WA(IW))
         GO TO 114
  113    CALL VRADBG (M,IDO,IP,L1,IDL1,CH,CH,CH,C,C,MDIMC,WA(IW))
  114    IF (IDO .EQ. 1) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDO
  116 CONTINUE
      SCALE=SQRT(1./N)
      IF (NA .EQ. 0) GO TO 118
      DO 117 J=1,N
      DO 117 I=1,M
         C(I,J) = SCALE*CH(I,J)
  117 CONTINUE
      RETURN
  118 DO 119 J=1,N
      DO 119 I=1,M
         C(I,J)=SCALE*C(I,J)
  119 CONTINUE
      RETURN
      END

!     ----------------------------------------------------------------- 
      SUBROUTINE VSCOSB(F,L,M,N,FT,C1,C2,WORK)
C
C     VSFFTPK, VERSION 1, AUGUST 1985
C
      integer i,j,k,k1,k2,l,m,n
      real F(L,M*N),FT(L,N,M),C1(M),C2(M),WORK(M+15)
      real scale

!     type declarations by j.w.lavelle, december, 1995

C
C     PREPROCESSING
C
      DO 100 K=1,N
         K1=M*(K-1)+1
         DO 100 I=1,L
  100       FT(I,K,1)=.5*F(I,K1)
      DO 200 J=2,M
      DO 200 K=1,N
         K1=M*(K-1)+J
         K2=M*K+2-J
         DO 200 I=1,L
  200       FT(I,K,J) =(C1(J)*F(I,K1)+C2(J)*F(I,K2))
C
C     REAL,PERIODIC ANALYSIS
C
      CALL VRFFTF(L*N,M,FT,F,L*N,WORK)
C
C     POSTPROCESSING
C
      SCALE=SQRT(2.)
      DO 320 K=1,N
         K1=N*(M-1)+K
         DO 300 I=1,L
  300       F(I,K) = SCALE*FT(I,K,1)
      IF (2*(M/2) .NE. M) GO TO 320
      DO 310 I=1,L
310         F(I,K1) = SCALE*FT(I,K,M)
  320 CONTINUE
      DO 400 J=2,M-1,2
      DO 400 K=1,N
         K1=N*(J-1)+K
         K2=K1+N
         DO 400 I=1,L
            F(I,K1)=SCALE*(FT(I,K,J)-FT(I,K,J+1))
400         F(I,K2)=SCALE*(FT(I,K,J)+FT(I,K,J+1))
      RETURN
      END

!     --------------------------------------------------------------------------

      SUBROUTINE VSCOSQ(F,L,M,N,FT,C1,C2,C3,C4,WORK)
C
C     VSFFTPK, VERSION 1, AUGUST 1985
C
      integer i,j,jby2,k,k1,k2,l,m,n
      real F(L,M*N),FT(L,N,M),C1(M),C2(M),C3(M),C4(M),WORK(M+15)

!     type declarations by j.w.lavelle, december, 1995
C
C     PREPROCESSING
C
      DO 120 K=1,N
         K1=M*(K-1)+1
         K2=M*K
         DO 100 I=1,L
  100       FT(I,K,1)=F(I,K1)
      IF (2*(M/2) .NE. M) GO TO 120
      DO 110 I=1,L
  110       FT(I,K,M)=-F(I,K2)
  120 CONTINUE
      DO 200 J=2,M-1,2
      JBY2=J/2
      DO 200 K=1,N
         K1=M*(K-1)+J
         K2=K1+1
         DO 200 I=1,L
            FT(I,K,J)   = F(I,K2)*C1(JBY2)+F(I,K1)*C2(JBY2)
  200       FT(I,K,J+1) = -F(I,K2)*C2(JBY2)+F(I,K1)*C1(JBY2)
C
C     REAL,PERIODIC SYNTHESIS
C
      CALL VRFFTB(L*N,M,FT,F,L*N,WORK)
C
C     POSTPROCESSING
C
      DO 300 J=1,M
      DO 300 K=1,N
         K1=N*(J-1)+K
         DO 300 I=1,L
  300       F(I,K1)=C4(J)*FT(I,K,J)+C3(J)*FT(I,K,M+1-J)
      RETURN
      END

!     ------------------------------------------------------------------------------

      SUBROUTINE VSRFTF(F,L,M,N,FT,WSAVE)
C
C     VSFFTPK, VERSION 1, AUGUST 1985
C
      integer i,j,k,k1,k2,l,m,n
      real F(L,M*N),FT(L,N,M),WSAVE(M+15)

!     type declarations by j.w.lavelle, december, 1995

C
C     RE-ORDER INPUT
C
      DO 100 K=1,N
      DO 100 J=1,M
      K1=M*(K-1)+J
      DO 100 I=1,L
  100 FT(I,K,J)=F(I,K1)
C
C     REAL, PERIODIC TRANSFORM
C
      CALL VRFFTF(L*N,M,FT,F,L*N,WSAVE)
      DO 220 K=1,N
      K1=N*(M-1)+K
      DO 200 I=1,L
  200 F(I,K)=FT(I,K,1)
      IF (2*(M/2) .NE. M) GO TO 220
      DO 210 I=1,L
  210 F(I,K1)=FT(I,K,M)
  220 CONTINUE
      DO 300 K=1,N
      DO 300 J=2,M-1,2
      K1=N*(J-1)+K
      K2=K1+N
      DO 300 I=1,L
      F(I,K1)=FT(I,K,J)
  300 F(I,K2)=FT(I,K,J+1)
      RETURN
      END

!     --------------------------------------------------------

      SUBROUTINE VSSINF(F,L,M,N,FT,C1,C2,WORK)
C
C     VSFFTPK, VERSION 1, AUGUST 1985
C
      integer i,j,k,k1,k2,l,m,n
      real F(L,M*N),FT(L,N,M),C1(M),C2(M),WORK(M+15)
      real scale

!     type declarations by j.w.lavelle, december, 1995
C
C     PREPROCESSING
C
      SCALE=SQRT(2.)
      DO 120 K=1,N
         K1=M*(K-1)+1
         K2=M*K
         DO 100 I=1,L
  100       FT(I,K,1)=SCALE*F(I,K1)
      IF (2*(M/2) .NE. M) GO TO 120
      DO 110 I=1,L
  110       FT(I,K,M)=-SCALE*F(I,K2)
  120 CONTINUE
      SCALE=.5*SCALE
      DO 200 J=2,M-1,2
      DO 200 K=1,N
         K1=M*(K-1)+J
         K2=K1+1
         DO 200 I=1,L
            FT(I,K,J)   = SCALE*(F(I,K2)-F(I,K1))
  200       FT(I,K,J+1) =-SCALE*(F(I,K2)+F(I,K1))
C
C     REAL,PERIODIC SYNTHESIS
C
      CALL VRFFTB(L*N,M,FT,F,L*N,WORK)
C
C     POSTPROCESSING
C
      DO 300 J=2,M
      DO 300 K=1,N
         K1=N*(J-2)+K
         DO 300 I=1,L
  300       F(I,K1)=C1(J)*FT(I,K,J)+C2(J)*FT(I,K,M+2-J)
      DO 400 K=1,N
         K1=N*(M-1)+K
         DO 400 I=1,L
  400       F(I,K1)=FT(I,K,1)
      RETURN
      END

!     ------------------------------------------------------------------

      SUBROUTINE VRADB4 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2,WA3)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      integer i,ic,ido,idp2,k,l1,m,mdimc,mp
      real  CC(MDIMC,IDO,4,L1)  ,CH(MDIMC,IDO,L1,4)    ,
     1                WA1(IDO)  ,WA2(IDO)  ,WA3(IDO)
      real sqrt2

!     type declarations by j.w.lavelle, december, 1995

      SQRT2=SQRT(2.)
      DO 101 K=1,L1
          DO 1001 M=1,MP
         CH(M,1,K,3) = (CC(M,1,1,K)+CC(M,IDO,4,K))
     1   -(CC(M,IDO,2,K)+CC(M,IDO,2,K))
         CH(M,1,K,1) = (CC(M,1,1,K)+CC(M,IDO,4,K))
     1   +(CC(M,IDO,2,K)+CC(M,IDO,2,K))
         CH(M,1,K,4) = (CC(M,1,1,K)-CC(M,IDO,4,K))
     1   +(CC(M,1,3,K)+CC(M,1,3,K))
         CH(M,1,K,2) = (CC(M,1,1,K)-CC(M,IDO,4,K))
     1   -(CC(M,1,3,K)+CC(M,1,3,K))
 1001     CONTINUE
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
               DO 1002 M=1,MP
            CH(M,I-1,K,1) = (CC(M,I-1,1,K)+CC(M,IC-1,4,K))
     1      +(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
            CH(M,I,K,1) = (CC(M,I,1,K)-CC(M,IC,4,K))
     1      +(CC(M,I,3,K)-CC(M,IC,2,K))
            CH(M,I-1,K,2)=WA1(I-2)*((CC(M,I-1,1,K)-CC(M,IC-1,4,K))
     1      -(CC(M,I,3,K)+CC(M,IC,2,K)))-WA1(I-1)
     1      *((CC(M,I,1,K)+CC(M,IC,4,K))+(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))
            CH(M,I,K,2)=WA1(I-2)*((CC(M,I,1,K)+CC(M,IC,4,K))
     1      +(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))+WA1(I-1)
     1      *((CC(M,I-1,1,K)-CC(M,IC-1,4,K))-(CC(M,I,3,K)+CC(M,IC,2,K)))
            CH(M,I-1,K,3)=WA2(I-2)*((CC(M,I-1,1,K)+CC(M,IC-1,4,K))
     1      -(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))-WA2(I-1)
     1      *((CC(M,I,1,K)-CC(M,IC,4,K))-(CC(M,I,3,K)-CC(M,IC,2,K)))
            CH(M,I,K,3)=WA2(I-2)*((CC(M,I,1,K)-CC(M,IC,4,K))
     1      -(CC(M,I,3,K)-CC(M,IC,2,K)))+WA2(I-1)
     1      *((CC(M,I-1,1,K)+CC(M,IC-1,4,K))-(CC(M,I-1,3,K)
     1      +CC(M,IC-1,2,K)))
            CH(M,I-1,K,4)=WA3(I-2)*((CC(M,I-1,1,K)-CC(M,IC-1,4,K))
     1      +(CC(M,I,3,K)+CC(M,IC,2,K)))-WA3(I-1)
     1     *((CC(M,I,1,K)+CC(M,IC,4,K))-(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))
            CH(M,I,K,4)=WA3(I-2)*((CC(M,I,1,K)+CC(M,IC,4,K))
     1      -(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))+WA3(I-1)
     1      *((CC(M,I-1,1,K)-CC(M,IC-1,4,K))+(CC(M,I,3,K)+CC(M,IC,2,K)))
 1002          CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
               DO 1003 M=1,MP
         CH(M,IDO,K,1) = (CC(M,IDO,1,K)+CC(M,IDO,3,K))
     1   +(CC(M,IDO,1,K)+CC(M,IDO,3,K))
         CH(M,IDO,K,2) = SQRT2*((CC(M,IDO,1,K)-CC(M,IDO,3,K))
     1   -(CC(M,1,2,K)+CC(M,1,4,K)))
         CH(M,IDO,K,3) = (CC(M,1,4,K)-CC(M,1,2,K))
     1   +(CC(M,1,4,K)-CC(M,1,2,K))
         CH(M,IDO,K,4) = -SQRT2*((CC(M,IDO,1,K)-CC(M,IDO,3,K))
     1   +(CC(M,1,2,K)+CC(M,1,4,K)))
 1003          CONTINUE
  106 CONTINUE
  107 RETURN
      END

!     -------------------------------------------------------------------------------

      SUBROUTINE VRADF2 (MP,IDO,L1,CC,CH,MDIMC,WA1)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      integer i,ic,ido,idp2,k,l1,m,mdimc,mp
      real   CH(MDIMC,IDO,2,L1)  ,CC(MDIMC,IDO,L1,2)     ,
     1                WA1(IDO)

!     type declarations by j.w.lavelle, december, 1995

      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = CC(M,1,K,1)+CC(M,1,K,2)
         CH(M,IDO,2,K) = CC(M,1,K,1)-CC(M,1,K,2)
 1001    CONTINUE
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            DO 1003 M=1,MP
            CH(M,I,1,K) = CC(M,I,K,1)+(WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))
            CH(M,IC,2,K) = (WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))-CC(M,I,K,1)
            CH(M,I-1,1,K) = CC(M,I-1,K,1)+(WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))
            CH(M,IC-1,2,K) = CC(M,I-1,K,1)-(WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))
 1003       CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         DO 1006 M=1,MP
         CH(M,1,2,K) = -CC(M,IDO,K,2)
         CH(M,IDO,1,K) = CC(M,IDO,K,1)
 1006    CONTINUE
  106 CONTINUE
  107 RETURN
      END

!     -------------------------------------------------------------

      SUBROUTINE VRADF5 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2,WA3,WA4)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      integer i,ic,ido,idp2,k,l1,m,mdimc,mp      
      real  CC(MDIMC,IDO,L1,5)    ,CH(MDIMC,IDO,5,L1)     ,
     1           WA1(IDO)     ,WA2(IDO)     ,WA3(IDO)     ,WA4(IDO)
      real arg,ti11,ti12,tr11,tr12,pimach

!     type declarations by j.w.lavelle, december, 1995

      ARG=2.*PIMACH(1.0)/5.
      TR11=COS(ARG)
      TI11=SIN(ARG)
      TR12=COS(2.*ARG)
      TI12=SIN(2.*ARG)
      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = CC(M,1,K,1)+(CC(M,1,K,5)+CC(M,1,K,2))+
     1    (CC(M,1,K,4)+CC(M,1,K,3))
         CH(M,IDO,2,K) = CC(M,1,K,1)+TR11*(CC(M,1,K,5)+CC(M,1,K,2))+
     1    TR12*(CC(M,1,K,4)+CC(M,1,K,3))
         CH(M,1,3,K) = TI11*(CC(M,1,K,5)-CC(M,1,K,2))+TI12*
     1    (CC(M,1,K,4)-CC(M,1,K,3))
         CH(M,IDO,4,K) = CC(M,1,K,1)+TR12*(CC(M,1,K,5)+CC(M,1,K,2))+
     1    TR11*(CC(M,1,K,4)+CC(M,1,K,3))
         CH(M,1,5,K) = TI12*(CC(M,1,K,5)-CC(M,1,K,2))-TI11*
     1    (CC(M,1,K,4)-CC(M,1,K,3))
 1001    CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DO 1002 M=1,MP
            CH(M,I-1,1,K) = CC(M,I-1,K,1)+((WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*
     1       CC(M,I,K,5)))+((WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))+(WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4)))
            CH(M,I,1,K) = CC(M,I,K,1)+((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))
            CH(M,I-1,3,K) = CC(M,I-1,K,1)+TR11*
     1      ( WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2)
     1       +WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5))+TR12*
     1      ( WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3)
     1       +WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))+TI11*
     1      ( WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2)
     1       -(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*CC(M,I-1,K,5)))+TI12*
     1      ( WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3)
     1       -(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*CC(M,I-1,K,4)))
            CH(M,IC-1,2,K) = CC(M,I-1,K,1)+TR11*
     1      ( WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2)
     1       +WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5))+TR12*
     1     ( WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3)
     1      +WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))-(TI11*
     1      ( WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2)
     1       -(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*CC(M,I-1,K,5)))+TI12*
     1      ( WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3)
     1       -(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*CC(M,I-1,K,4))))
            CH(M,I,3,K) = (CC(M,I,K,1)+TR11*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR12*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))+(TI11*((WA4(I-2)*CC(M,I-1,K,5)+
     1       WA4(I-1)*CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))+TI12*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))
            CH(M,IC,2,K) = (TI11*((WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*
     1       CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))+TI12*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))-(CC(M,I,K,1)+TR11*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR12*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
            CH(M,I-1,5,K) = (CC(M,I-1,K,1)+TR12*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*
     1       CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5)))+TR11*((WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))+(WA3(I-2)*
     1       CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))))+(TI12*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA4(I-2)*CC(M,I,K,5)-
     1       WA4(I-1)*CC(M,I-1,K,5)))-TI11*((WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
            CH(M,IC-1,4,K) = (CC(M,I-1,K,1)+TR12*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*
     1       CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5)))+TR11*((WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))+(WA3(I-2)*
     1       CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))))-(TI12*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA4(I-2)*CC(M,I,K,5)-
     1       WA4(I-1)*CC(M,I-1,K,5)))-TI11*((WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
            CH(M,I,5,K) = (CC(M,I,K,1)+TR12*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR11*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))+(TI12*((WA4(I-2)*CC(M,I-1,K,5)+
     1       WA4(I-1)*CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))-TI11*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))
            CH(M,IC,4,K) = (TI12*((WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*
     1       CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))-TI11*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))-(CC(M,I,K,1)+TR12*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR11*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
 1002       CONTINUE
  102    CONTINUE
  103 CONTINUE
      RETURN
      END

!     -------------------------------------------------------------

      SUBROUTINE VRFFTF (M,N,R,RT,MDIMR,WSAVE)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      integer m,mdimr,n
      real       R(MDIMR,N)  ,RT(MDIMR,N)    ,WSAVE(N+15)

!     type declarations by j.w.lavelle, december, 1995

      IF (N .EQ. 1) RETURN
      CALL VRFTF1 (M,N,R,RT,MDIMR,WSAVE(1),WSAVE(N+1))
      RETURN
      END

!     ------------------------------------------------------------

      SUBROUTINE VRFTF1 (M,N,C,CH,MDIMC,WA,FAC)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      integer i,idl1,ido,ip,iw,ix2,ix3,ix4,j,k1,kh,l1,l2,m,mdimc,n,na,nf
      real       CH(MDIMC,N) ,C(MDIMC,N)  ,WA(N)   ,FAC(15)
      real scale
!     type declarations by j.w.lavelle, december, 1995

      NF = FAC(2)
      NA = 1
      L2 = N
      IW = N
      DO 111 K1=1,NF
         KH = NF-K1
         IP = FAC(KH+3)
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW-(IP-1)*IDO
         NA = 1-NA
         IF (IP .NE. 4) GO TO 102
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL VRADF4 (M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  101    CALL VRADF4 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  102    IF (IP .NE. 2) GO TO 104
         IF (NA .NE. 0) GO TO 103
         CALL VRADF2 (M,IDO,L1,C,CH,MDIMC,WA(IW))
         GO TO 110
  103    CALL VRADF2 (M,IDO,L1,CH,C,MDIMC,WA(IW))
         GO TO 110
  104    IF (IP .NE. 3) GO TO 106
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 105
         CALL VRADF3 (M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2))
         GO TO 110
  105    CALL VRADF3 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2))
         GO TO 110
  106    IF (IP .NE. 5) GO TO 108
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 107
      CALL VRADF5(M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  107 CALL VRADF5 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  108    IF (IDO .EQ. 1) NA = 1-NA
         IF (NA .NE. 0) GO TO 109
         CALL VRADFG (M,IDO,IP,L1,IDL1,C,C,C,CH,CH,MDIMC,WA(IW))
         NA = 1
         GO TO 110
  109    CALL VRADFG (M,IDO,IP,L1,IDL1,CH,CH,CH,C,C,MDIMC,WA(IW))
         NA = 0
  110    L2 = L1
  111 CONTINUE
      SCALE=SQRT(1./N)
      IF (NA .EQ. 1) GO TO 113
      DO 112 J=1,N
      DO 112 I=1,M
         C(I,J) = SCALE*CH(I,J)
  112 CONTINUE
      RETURN
  113 DO 114 J=1,N
      DO 114 I=1,M
         C(I,J)=SCALE*C(I,J)
  114 CONTINUE
      RETURN
      END

!     ---------------------------------------------------------------

      SUBROUTINE VSCOSF(F,L,M,N,FT,C1,C2,WORK)
C
C     VSFFTPK, VERSION 1, AUGUST 1985
C
      integer i,j,k,k1,k2,l,m,n
      real F(L,N*M),FT(L,N,M),C1(M),C2(M),WORK(M+15)
      real scale

!     type declarations by j.w.lavelle, december, 1995

C
C     PREPROCESSING
C
      SCALE=SQRT(2.)
      DO 120 K=1,N
         K1=M*(K-1)+1
         K2=M*K
         DO 100 I=1,L
  100       FT(I,K,1)=SCALE*F(I,K1)
      IF (2*(M/2) .NE. M) GO TO 120
      DO 110 I=1,L
  110       FT(I,K,M)=SCALE*F(I,K2)
  120 CONTINUE
      SCALE=.5*SCALE
      DO 200 J=2,M-1,2
      DO 200 K=1,N
         K1=M*(K-1)+J
         K2=K1+1
         DO 200 I=1,L
            FT(I,K,J)   = SCALE*(F(I,K1)+F(I,K2))
  200       FT(I,K,J+1) =-SCALE*(F(I,K1)-F(I,K2))
C
C     REAL,PERIODIC SYNTHESIS
C
      CALL VRFFTB(L*N,M,FT,F,L*N,WORK)
C
C     POSTPROCESSING
C
      DO 300 J=2,M
      DO 300 K=1,N
         K1=N*(J-1)+K
         DO 300 I=1,L
  300       F(I,K1)=C1(J)*FT(I,K,J)-C2(J)*FT(I,K,M+2-J)
      DO 400 K=1,N
         DO 400 I=1,L
  400       F(I,K)=FT(I,K,1)
      RETURN
      END

!     --------------------------------------------------------------------

      SUBROUTINE VSCSQI(N,C1,C2,C3,C4,WSAVE)
C
C     VSFFTPK, VERSION 1, AUGUST 1985
C
      integer i,n
      ENTRY VSSNQI(N,C1,C2,C3,C4,WSAVE)
      real C1(N),C2(N),C3(N),C4(N),WSAVE(N+15)
      real c,dx,pi,s,scale,pimach

!     type declarations by j.w.lavelle, december, 1995

      PI=PIMACH(1.0)
      DX=PI/N
      SCALE=SQRT(.5)
C
C     GENERATE A(I)+-B(I)
C
      DO 100 I=1,(N-1)/2
         C=COS(I*DX)
         S=SIN(I*DX)
         C1(I)=.5*(S+C)
  100    C2(I)=.5*(C-S)
C
      DX=PI/(2*N)
      DO 200 I=1,N
         C=COS((I-.5)*DX)
         S=SIN((I-.5)*DX)
         C3(I)=SCALE*(C+S)
200      C4(I)=SCALE*(C-S)
C
C     INITIALIZE VRFFTPK ROUTINE
C
      CALL VRFFTI(N,WSAVE)
      RETURN
      END
!       end  SUBROUTINE VSCSQI
!     --------------------------------------------------------------

      SUBROUTINE VSRFTI(N,WSAVE)
C
C     VSFFTPK, VERSION 1, AUGUST 1985
C
      integer n
      real WSAVE(N+15)

!     type declarations by j.w.lavelle, december, 1995
C
C     INITIALIZE VRFFTPK ROUTINE
C
      CALL VRFFTI(N,WSAVE)
      RETURN
      END

!     ------------------------------------------------------------------------

!      FUNCTION PIMACH (DUM)
       real function pimach(dum)
C***BEGIN PROLOGUE  PIMACH
C***SUBSIDIARY
C***PURPOSE  Subsidiary to HSTCSP, HSTSSP and HWSCSP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (PIMACH-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subprogram supplies the value of the constant PI correct to
C     machine precision where
C
C     PI=3.1415926535897932384626433832795028841971693993751058209749446
C
C***SEE ALSO  HSTCSP, HSTSSP, HWSCSP
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  PIMACH
C
C***FIRST EXECUTABLE STATEMENT  PIMACH

!     real pimach,dum
      real  dum
!     type declarations by j.w.lavelle, december, 1995

      PIMACH = 3.1415926535897932
      RETURN
      END

!     -----------------------------------------------------
      SUBROUTINE VRFFTI (N,WSAVE)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      integer n
      real       WSAVE(N+15)

!     type declarations by j.w.lavelle, december, 1995

      IF (N .EQ. 1) RETURN
      CALL VRFTI1 (N,WSAVE(1),WSAVE(N+1))
      RETURN
      END

!     --------------------------------------------------------------------------

      SUBROUTINE VRFTI1 (N,WA,FAC)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      integer i,ib,ido,ii,ip,ipm,is,j,k1,l1,l2,ld,n,nf,nfm1,nl,nq,nr,ntry,NTRYH(4)
      real WA(N),FAC(15)  
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
      real arg,argh,argld,fi,tpi,pimach 

!     type declarations by j.w.lavelle, december, 1995

      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      FAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         FAC(IB+2) = FAC(IB+1)
  106 CONTINUE
      FAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      FAC(1) = N
      FAC(2) = NF
      TPI = 2.*PIMACH(1.0)
      ARGH = TPI/REAL(N)
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
      DO 110 K1=1,NFM1
         IP = FAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         DO 109 J=1,IPM
            LD = LD+L1
            I = IS
            ARGLD = REAL(LD)*ARGH
            FI = 0.
            DO 108 II=3,IDO,2
               I = I+2
               FI = FI+1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IS = IS+IDO
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END
