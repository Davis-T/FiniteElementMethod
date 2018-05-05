PROGRAM Finite_element_Analysis
    IMPLICIT NONE
    INTEGER::I
    INTEGER::NN,NE,ND,NFIX,NFORCE,NTYPE !NN结点总数,NE单元总数,ND总自由度数,NFIX约束数,NFORCE荷载总数,NTYPE求解类型（1平面应力问题、2平面应变问题）
    REAL::E,ANU,T,GM                    !E弹性模量,ANU泊松比,T结构厚度,GM材料重度
    REAL*8,DIMENSION(:),ALLOCATABLE::CX,CY,LFIX,FN,F,DISP     !CX,CY结点坐标,LFIX约束编号,FN荷载信息,F结点荷载列阵,DISP结点位移列阵
    REAL*8,DIMENSION(:,:),ALLOCATABLE::LOC,TK              !LOC单元结点编号,TK总刚
    OPEN(5,FILE='DATA.txt',STATUS='OLD')
    OPEN(6,FILE='OUT.txt',STATUS='NEW')
    OPEN(7,FILE='TK.txt',STATUS='NEW')
    OPEN(8,FILE='AREA.txt',STATUS='NEW')
    OPEN(9,FILE='B.txt',STATUS='NEW')
    OPEN(10,FILE='D.txt',STATUS='NEW')
    OPEN(11,FILE='KE.txt',STATUS='NEW')
    READ(5,*)NN,NE,NFIX,E,ANU,T,GM,NFORCE,NTYPE     !读取控制信息
    WRITE(6,*)NN,NE,NFIX,E,ANU,T,GM,NFORCE,NTYPE
    ND=NN*2
    WRITE(6,*)"ND=",ND
    IF(NTYPE /= 1) THEN
        E=E/(1.0-ANU*ANU)
        ANU=ANU/(1.0-ANU)
    END IF
    ALLOCATE(CX(NN),CY(NN),LFIX(NFIX),LOC(NE,3),TK(ND,ND),FN(NFORCE*2),F(ND),DISP(ND))
    CALL DATAIN(NN,NE,NFIX,NFORCE,LOC,CX,CY,LFIX,FN)   !读取结构信息——LOC单元结点编码、CX,CY结点坐标、LFIX约束信息、FN荷载信息
    CALL GLOSTIF(NN,NE,ND,E,ANU,T,LOC,CX,CY,TK)    !计算总刚矩阵TK
    CALL GLOFORC(NN,NE,ND,NFORCE,T,GM,LOC,CX,CY,FN,F)    !计算结点荷载列阵F
    CALL INSCD(ND,NFIX,LFIX,TK,F,DISP)     !引入边界条件
    CALL SOLVE(ND,TK,DISP)     !求解位移DISP
    CALL INNERFORCE(NN,NE,ND,E,ANU,T,GM,LOC,CX,CY,DISP)    !计算内力
    CLOSE(5)
    CLOSE(6)
END PROGRAM Finite_element_Analysis

!==================输入数据==================
SUBROUTINE DATAIN(NN,NE,NFIX,NFORCE,LOC,CX,CY,LFIX,FN)  
!NN结点总数,NE单元总数,NFIX约束数,NFORCE荷载总数,LOC单元编码,CX,CY编码坐标,LFIX约束信息,FN荷载信息
    IMPLICIT NONE
    INTEGER::I,J,K,L
    INTEGER::NN,NE,NFIX,NFORCE
    REAL*8,DIMENSION(:)::CX(NN),CY(NN),LFIX(NFIX),FN(NFORCE)
    REAL*8,DIMENSION(:,:)::LOC(NE,3)
    READ(5,*)(LOC(I,1),LOC(I,2),LOC(I,3),I=1,NE)
    READ(5,*)(CX(J),CY(J),J=1,NN)
    READ(5,*)(LFIX(K),K=1,NFIX)
    READ(5,*)(FN(L),L=1,NFORCE*2)
END
!==================输入数据==================

!================计算单元面积=================
SUBROUTINE AREA(IE,NN,NE,LOC,CX,CY,AE)
!IE单元计数,NN结点总数,NE单元总数,LOC结点编码,CX,CY结点坐标,AE单元面积
    IMPLICIT NONE
    INTEGER::I,J,K
    INTEGER::IE,NN,NE
    REAL::AE,XJI,XKI,YJI,YKI
    REAL*8,DIMENSION(:)::CX(NN),CY(NN)
    REAL*8,DIMENSION(:,:)::LOC(NE,3)
    I=LOC(IE,1)
    J=LOC(IE,2)
    K=LOC(IE,3)
    XJI=CX(J)-CX(I)
    XKI=CX(K)-CX(I)
    YJI=CY(J)-CY(I)
    YKI=CY(K)-CY(I)
    AE=0.5*(XJI*YKI-XKI*YJI)
    WRITE(8,'(I3,F8.4)')IE,AE
END
!================计算单元面积================

!==============计算单元应变矩阵===============
SUBROUTINE CALMB(IE,NN,NE,LOC,CX,CY,AE,MB)
!IE单元计数,NN结点总数,NE单元总数,LOC单元结点编码,CX,CY结点坐标,AE单元面积,MB应变矩阵
    IMPLICIT NONE
    INTEGER::I,J,K,II,JJ,KK     !II,JJ,KK循环计数
    INTEGER::IE,NN,NE
    REAL::AE
    REAL*8,DIMENSION(:)::CX(NN),CY(NN)
    REAL*8,DIMENSION(:,:)::LOC(NE,3),MB(3,6)
    CALL AREA(IE,NN,NE,LOC,CX,CY,AE)    !计算单元面积AE
    I=LOC(IE,1)
    J=LOC(IE,2)
    K=LOC(IE,3)
    DO II=1,3
        DO JJ=1,6
            MB(II,JJ)=0.0
        END DO
    END DO
    MB(1,1)=CY(J)-CY(K)
    MB(1,3)=CY(K)-CY(I)
    MB(1,5)=CY(I)-CY(J)
    MB(2,2)=-CX(J)+CX(K)
    MB(2,4)=-CX(K)+CX(I)
    MB(2,6)=-CX(I)+CX(J)
    MB(3,1)=MB(2,2)
    MB(3,2)=MB(1,1)
    MB(3,3)=MB(2,4)
    MB(3,4)=MB(1,3)
    MB(3,5)=MB(2,6)
    MB(3,6)=MB(1,5)
    DO II=1,3
        DO JJ=1,6
            MB(II,JJ)=MB(II,JJ)/(2*AE)
        END DO
    END DO
    DO I=1,3
        WRITE(9,'(I3,6F8.4)')IE,(MB(I,J),J=1,6)
    END DO
END
!===============计算单元应变矩阵===============

!================计算弹性矩阵=================
SUBROUTINE CALMD(E,ANU,MD)
!E弹性模量,ANU泊松比,MD弹性矩阵
    IMPLICIT NONE
    INTEGER::I,J
    REAL::E,ANU
    REAL*8,DIMENSION(:,:)::MD(3,3)
    DO I=1,3
        DO J=1,3
            MD(I,J)=0.0
        END DO
    END DO
    MD(1,1)=E/(1.0-ANU*ANU)
    MD(1,2)=ANU*E/(1.0-ANU*ANU)
    MD(2,2)=MD(1,1)
    MD(2,1)=MD(1,2)
    MD(3,3)=0.5*E/(1.0+ANU)
    DO I=1,3
        WRITE(10,'(3F8.4)')(MD(I,J),J=1,3)
    END DO
END
!================计算弹性矩阵================

!===============计算单元刚度矩阵===============
SUBROUTINE STE(IE,NN,NE,E,ANU,T,LOC,CX,CY,KE)    
!IE单元计数,NN结点总数,NE单元总数,E弹性模量,ANU泊松比,T单元厚度,LOC结点编码,CX,CY结点坐标,KE单元刚度矩阵
    IMPLICIT NONE
    INTEGER::IE,NN,NE,I,J
    REAL::E,ANU,T,AE=0.0
    REAL*8,DIMENSION(:)::CX(NN),CY(NN)
    REAL*8,DIMENSION(:,:)::LOC(NE,3),KE(6,6),MB(3,6),MD(3,3),MBT(6,3),MS(3,6)   !MB应变矩阵,MD弹性矩阵,MBT应变矩阵的转置,MS应力矩阵
    DO I=1,3
        DO J=1,6
            MB(I,J)=0.0
        END DO
    END DO
    DO I=1,3
        DO J=1,3
            MD(I,J)=0.0
        END DO
    END DO
    CALL AREA(IE,NN,NE,LOC,CX,CY,AE)
    CALL CALMB(IE,NN,NE,LOC,CX,CY,AE,MB)
    CALL CALMD(E,ANU,MD)
    MBT=TRANSPOSE(MB)
    MS=MATMUL(MD,MB)
    KE=MATMUL(MBT,MS)
    DO I=1,6
        DO J=1,6
            KE(I,J)=KE(I,J)*T*AE
        END DO
    END DO
    DO I=1,6
        WRITE(11,'(I3,6F8.4)')IE,(KE(I,J),J=1,6)
    END DO
END
!===============计算单元刚度矩阵===============

!===============集成总体刚度矩阵===============
SUBROUTINE GLOSTIF(NN,NE,ND,E,ANU,T,LOC,CX,CY,TK)
!NN结点总数,NE单元总数,ND总自由度,E弹性模量,ANU泊松比,T厚度,LOC单元结点编码,CX,CY节点坐标,TK总刚
    IMPLICIT NONE
    INTEGER::I,J,II,JJ,IZ,JZ    !IZ,JZ定位编码
    INTEGER::IE,NN,NE,ND    !IE单元计数,NN结点总数,NE单元总数,ND总自由度
    REAL::E,ANU,T   !E弹性模量,ANU泊松比,T单元结构厚度
    REAL*8,DIMENSION(:)::CX(NN),CY(NN)  !结点坐标
    REAL*8,DIMENSION(:,:)::KE(6,6),LOC(NE,3),TK(ND,ND)  !KE单刚,LOC结点编码,TK总刚
    DO I=1,ND
        DO J=1,ND
            TK(I,J)=0.0
        END DO
    END DO
    DO IE=1,NE
        CALL STE(IE,NN,NE,E,ANU,T,LOC,CX,CY,KE)
        DO I=1,6
            IF(I==1.OR.I==2) THEN
                II=LOC(IE,1)
            ELSE IF(I==3.OR.I==4) THEN
                II=LOC(IE,2)
            ELSE
                II=LOC(IE,3)
            ENDIF
            DO J=1,6
                IF(J==1.OR.J==2) THEN
                    JJ=LOC(IE,1)
                ELSE IF(J==3.OR.J==4) THEN
                    JJ=LOC(IE,2)
                ELSE
                    JJ=LOC(IE,3)
                ENDIF
                IF(MOD(I,2)==0) THEN
                    IZ=II*2
                ELSE
                    IZ=II*2-1
                ENDIF
                IF(MOD(J,2)==0) THEN
                    JZ=JJ*2
                ELSE
                    JZ=JJ*2-1
                ENDIF
                TK(IZ,JZ)=TK(IZ,JZ)+KE(I,J)
            END DO
        END DO
    END DO
END
!===============集成总体刚度矩阵===============

!=================计算结点荷载================
SUBROUTINE GLOFORC(NN,NE,ND,NFORCE,T,GM,LOC,CX,CY,FN,F)    
!NN结点总数,NE单元总数,ND自由度总数,NFORCE荷载总数,T厚度,GM重度,LOC单元结点编码,CX,CY结点坐标,FG等效重力荷载,FN荷载信息,F总结点荷载列阵
    INTEGER::NN,NE,ND,IE,NFORCE
    INTEGER::I,II   !II集成结点荷载定位编码
    REAL::T,GM,FG=0.0,A=0.0,B=0.0
    REAL*8,DIMENSION(:)::CX(NN),CY(NN),FN(NFORCE*2),F(ND)
    REAL*8,DIMENSION(:,:)::LOC(NE,3)
    DO I=1,ND
        F(I)=0.0
    END DO
    DO I=1,NFORCE
        F(INT(FN(2*I-1)))=FN(2*I)
    END DO
    DO IE=1,NE
        CALL AREA(IE,NN,NE,LOC,CX,CY,AE)
        FG=-GM*T*AE/3
        DO I=1,3
            II=LOC(IE,I)
            F(2*II)=F(2*II)+FG     !外力结点荷载+重力等效结点荷载
        END DO
    END DO
END
!=================计算结点荷载================

!=================引入边界条件================
SUBROUTINE INSCD(ND,NFIX,LFIX,TK,F,DISP)
!ND总自由度,NFIX总约束数,LFIX约束信息,TK总刚,F荷载列阵,DISP结点位移
    IMPLICIT NONE
    INTEGER::ND,NFIX,I,J,II
    REAL*8,DIMENSION(:)::LFIX(NFIX),F(ND),DISP(ND)
    REAL*8,DIMENSION(:,:)::TK(ND,ND)
    DO I=1,NFIX
        II=LFIX(I)
        TK(II,II)=1.0
        F(II)=0.0
        DO J=1,ND
            IF(J/=II) THEN
                TK(II,J)=0.0
                TK(J,II)=0.0
            ENDIF
        END DO
    END DO
    DO I=1,ND
        DISP(I)=F(I)
        DO J=1,ND
        WRITE(7,'(F6.4,",",$)') TK(I,J)
        END DO
        WRITE(7,'(";",$)')
    END DO
END
!=================引入边界条件================

!==============高斯消去法求解方程==============
SUBROUTINE SOLVE(N,A,B)
    IMPLICIT NONE
    INTEGER::N,I,J,I1,M
    REAL*8,DIMENSION(:,:)::A(N,N)
    REAL*8,DIMENSION(:)::B(N)
    DO I=1,N
        I1=I+1
        DO J=I1,N
            A(I,J)=A(I,J)/A(I,I)
        END DO
        B(I)=B(I)/A(I,I)
        A(I,I)=1.0
        DO J=I1,N
            DO M=I1,N
                A(J,M)=A(J,M)-A(J,I)*A(I,M)
            END DO
            B(J)=B(J)-A(J,I)*B(I)
        END DO
    END DO
    DO I=N-1,1,-1
        DO J=I+1,N
            B(I)=B(I)-A(I,J)*B(J)
        END DO
    END DO
    DO I=1,N,2
        WRITE(6,*) B(I),B(I+1)
    END DO
END
!==============高斯消去法求解方程==============

!=================单元应力===================
SUBROUTINE INNERFORCE(NN,NE,ND,E,ANU,T,GM,LOC,CX,CY,DISP)
!IE单元计数,NN结点总数,NE单元总数,ND总自由度,E弹性模量,ANU泊松比,T厚度,GM重度,LOC单元结点编码,CX,CY结点坐标,DISP结点位移
    IMPLICIT NONE
    INTEGER::I,II   !I单元结点位移定位码,II总体位移定位码
    INTEGER::IE,NN,NE,ND
    REAL::E,ANU,T,GM,AE
    REAL*8,DIMENSION(:)::CX(NN),CY(NN),DISP(ND),INNF(3),DISPE(6)      !DISP结点位移,INNF单元内力,DISPE单元结点位移
    REAL*8,DIMENSION(:,:)::LOC(NE,3),MB(3,6),MD(3,3),MS(3,6)      !MB应变矩阵,MD弹性矩阵,MS应力矩阵
    DO IE=1,NE
        CALL AREA(IE,NN,NE,LOC,CX,CY,AE)    !计算单元面积
        CALL CALMB(IE,NN,NE,LOC,CX,CY,AE,MB)    !计算单元应变矩阵B
        CALL CALMD(E,ANU,MD)    !计算弹性矩阵D
        MS=MATMUL(MD,MB)    !计算单元应力矩阵S
        DO I=1,6
            IF(MOD(I,2)==0) THEN
                II=LOC(IE,I/2)*2
            ELSE
                II=LOC(IE,(I+1)/2)*2-1
            ENDIF
            DISPE(I)=DISP(II)
        END DO
        INNF=MATMUL(MS,DISPE)
        WRITE(*,'(1X,I3,3E15.6)')IE,(INNF(I),I=1,3)
    END DO
END
!=================单元应力===================