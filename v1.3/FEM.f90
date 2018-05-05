PROGRAM Finite_element_Analysis
    IMPLICIT NONE
    INTEGER::I,J
    INTEGER::NN,NE,ND,NW,NFIX,NFORCE,NTYPE!NN结点总数,NE单元总数,ND总自由度数,NW半带宽,NFIX约束数,NFORCE荷载总数,NTYPE求解类型（1平面应力问题、2平面应变问题）
    REAL*8::E,ANU,T,GM,STARTTIME,ENDTIME,TOTLETIME           !E弹性模量,ANU泊松比,T结构厚度,GM材料重度
    REAL*8,DIMENSION(:),ALLOCATABLE::CX,CY,LFIX,FN,F,DISP     !CX,CY结点坐标,LFIX约束编号,FN荷载信息,F结点荷载列阵,DISP结点位移列阵
    REAL*8,DIMENSION(:,:),ALLOCATABLE::LOC,TK              !LOC单元结点编号,TK总刚
	CALL CPU_TIME(STARTTIME)
    OPEN(50,FILE='DATA.txt',STATUS='OLD')
    OPEN(60,FILE='INFO.txt',STATUS='NEW')
	OPEN(70,FILE='DISP.txt',STATUS='NEW')
	OPEN(80,FILE='STR.txt',STATUS='NEW')
    READ(50,*)NN,NE,NFIX,E,ANU,T,GM,NFORCE,NTYPE     !读取控制信息
    WRITE(60,100)
    WRITE(60,101)NN,NE,NFIX,E,ANU,T,GM,NFORCE,NTYPE
100 FORMAT(4X,"NN",4X,"NE",4X,"NFIX",7X,"E",8X,"ANU",6X,"T",6X,"GM",5X,"NFORCE",2X,"NTYPE") 
101 FORMAT(2X,I4,2X,I4,2X,I4,1X,E13.3,2X,F6.4,2X,F6.4,2X,F6.4,3X,I4,5X,I2)
    ND=NN*2
    IF(NTYPE /= 1) THEN
        E=E/(1.0-ANU*ANU)
        ANU=ANU/(1.0-ANU)
    END IF
    ALLOCATE(CX(NN),CY(NN),LFIX(NFIX),LOC(NE,3),FN(NFORCE*2),F(ND),DISP(ND))
    CALL DATAIN(NN,NE,NW,NFIX,NFORCE,LOC,CX,CY,LFIX,FN)   !读取结构信息——LOC单元结点编码、CX,CY结点坐标、LFIX约束信息、FN荷载信息
	ALLOCATE(TK(ND,NW))
    CALL GLOSTIF(NN,NE,ND,NW,E,ANU,T,LOC,CX,CY,TK)    !计算总刚矩阵TK
    CALL GLOFORC(NN,NE,ND,NFORCE,T,GM,LOC,CX,CY,FN,F)    !计算结点荷载列阵F
    CALL INSCD(ND,NW,NFIX,LFIX,TK,F,DISP)     !引入边界条件
    CALL SOLVE(ND,NW,TK,DISP)     !求解位移DISP
    WRITE(70,200)
    WRITE(70,201)(I,DISP(I*2-1),DISP(I*2),I=1,NN)
200 FORMAT(4X,"NODE",5X,"X-DISP",9X,"Y-DISP")
201 FORMAT(2X,I4,2X,E13.3,2X,E13.3)
    WRITE(70,'(/)')
    WRITE(80,300)
300 FORMAT(3X,"ELEMENT",7X,"X-STR",10X,"Y-STR",10X,"XY-STR")
    CALL INNERFORCE(NN,NE,ND,E,ANU,T,GM,LOC,CX,CY,DISP)    !计算内力并打印
	DEALLOCATE(CX,CY,LFIX,LOC,TK,FN,F,DISP)
    CALL CPU_TIME(ENDTIME)
	TOTLETIME=ENDTIME-STARTTIME
	WRITE(60,'(//,4X,"TotleTime:",F12.6,2X,"seconds")') TOTLETIME
	CLOSE(50)
    CLOSE(60)
	CLOSE(70)
	CLOSE(80)
	PRINT '(/,4X,A)', "AnalysisComplete!"
	READ(*,*)
END PROGRAM Finite_element_Analysis

!==================输入数据==================
SUBROUTINE DATAIN(NN,NE,NW,NFIX,NFORCE,LOC,CX,CY,LFIX,FN)  
!NN结点总数,NE单元总数,NFIX约束数,NFORCE荷载总数,LOC单元编码,CX,CY编码坐标,LFIX约束信息,FN荷载信息
    IMPLICIT NONE
	CHARACTER(LEN=50)::CTM="  1 ReadData:"
    INTEGER::I,J,K,L,IW,IE
    INTEGER::NN,NE,NW,NFIX,NFORCE
    REAL*8,DIMENSION(:)::CX(NN),CY(NN),LFIX(NFIX),FN(NFORCE)
    REAL*8,DIMENSION(:,:)::LOC(NE,3)
    READ(50,*)(LOC(I,1),LOC(I,2),LOC(I,3),I=1,NE)    !读取结点编码
    READ(50,*)(CX(J),CY(J),J=1,NN)       !读取结点坐标
    READ(50,*)(LFIX(K),K=1,NFIX)     !读取约束信息
    READ(50,*)(FN(L),L=1,NFORCE*2)       !读取荷载信息
	NW=0
    DO IE=1,NE
        DO I=1,2
            DO J=I+1,3
                IW=ABS(LOC(IE,I)-LOC(IE,J))
                IF(NW<IW)THEN
                    NW=IW
                END IF
            END DO
        END DO
    END DO
	NW=(NW+1)*2
	PRINT '(A,10X,"DONE!")',CTM
	RETURN
END
!==================输入数据==================

!================计算单元面积=================
SUBROUTINE AREA(IE,NN,NE,LOC,CX,CY,AE)
!IE单元计数,NN结点总数,NE单元总数,LOC结点编码,CX,CY结点坐标,AE单元面积
    IMPLICIT NONE
    INTEGER::I,J,K
    INTEGER::IE,NN,NE
    REAL*8::AE,XJI,XKI,YJI,YKI
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
	RETURN
END
!================计算单元面积================

!==============计算单元应变矩阵===============
SUBROUTINE CALMB(IE,NN,NE,LOC,CX,CY,AE,MB)
!IE单元计数,NN结点总数,NE单元总数,LOC单元结点编码,CX,CY结点坐标,AE单元面积,MB应变矩阵
    IMPLICIT NONE
    INTEGER::I,J,K,II,JJ,KK     !II,JJ,KK循环计数
    INTEGER::IE,NN,NE
    REAL*8::AE
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
	RETURN
END
!===============计算单元应变矩阵===============

!================计算弹性矩阵=================
SUBROUTINE CALMD(E,ANU,MD)
!E弹性模量,ANU泊松比,MD弹性矩阵
    IMPLICIT NONE
    INTEGER::I,J
    REAL*8::E,ANU
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
	RETURN
END
!================计算弹性矩阵================

!===============计算单元刚度矩阵===============
SUBROUTINE STE(IE,NN,NE,E,ANU,T,LOC,CX,CY,KE)    
!IE单元计数,NN结点总数,NE单元总数,E弹性模量,ANU泊松比,T单元厚度,LOC结点编码,CX,CY结点坐标,KE单元刚度矩阵
    IMPLICIT NONE
    INTEGER::IE,NN,NE,I,J
    REAL*8::E,ANU,T,AE=0.0
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
	RETURN
END
!===============计算单元刚度矩阵===============

!===============集成总体刚度矩阵===============
SUBROUTINE GLOSTIF(NN,NE,ND,NW,E,ANU,T,LOC,CX,CY,TK)
!NN结点总数,NE单元总数,ND总自由度,NW半带宽,E弹性模量,ANU泊松比,T厚度,LOC单元结点编码,CX,CY节点坐标,TK总刚
    IMPLICIT NONE
	CHARACTER(LEN=50)::CTM="  2 MakeGlobalStiffnessMatrix:"
    INTEGER::I,J,II,JJ    !II,JJ定位编码
    INTEGER::IE,NN,NE,ND,NW    !IE单元计数,NN结点总数,NE单元总数,ND总自由度
    REAL*8::E,ANU,T   !E弹性模量,ANU泊松比,T单元结构厚度
    REAL*8,DIMENSION(:)::CX(NN),CY(NN)  !结点坐标
    REAL*8,DIMENSION(:,:)::KE(6,6),LOC(NE,3),TK(ND,NW)  !KE单刚,LOC结点编码,TK总刚
    DO I=1,ND
        DO J=1,NW
            TK(I,J)=0.0
        END DO
    END DO
    DO IE=1,NE
        CALL STE(IE,NN,NE,E,ANU,T,LOC,CX,CY,KE)
        DO I=1,6
            IF(MOD(I,2)==0) THEN
                II=LOC(IE,I/2)*2
            ELSE
                II=LOC(IE,(I+1)/2)*2-1
            ENDIF
            DO J=1,6
                IF(MOD(J,2)==0) THEN
                    JJ=LOC(IE,J/2)*2
                ELSE
                    JJ=LOC(IE,(J+1)/2)*2-1
                ENDIF
				IF(JJ>=II) THEN
					JJ=JJ-II+1
					TK(II,JJ)=TK(II,JJ)+KE(I,J)
				END IF
            END DO
        END DO
    END DO
	PRINT '(A,10X,"DONE!")',CTM
	RETURN
END
!===============集成总体刚度矩阵===============

!=================计算结点荷载================
SUBROUTINE GLOFORC(NN,NE,ND,NFORCE,T,GM,LOC,CX,CY,FN,F)    
!NN结点总数,NE单元总数,ND自由度总数,NFORCE荷载总数,T厚度,GM重度,LOC单元结点编码,CX,CY结点坐标,FG等效重力荷载,FN荷载信息,F总结点荷载列阵
	IMPLICIT NONE
	CHARACTER(LEN=50)::CTM="  3 MakeLoadArray:"
	INTEGER::NN,NE,ND,IE,NFORCE
    INTEGER::I,II   !II集成结点荷载定位编码
    REAL*8::T,GM,FG=0.0,A=0.0,B=0.0,AE
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
	PRINT '(A,10X,"DONE!")',CTM
	RETURN
END
!=================计算结点荷载================

!=================引入边界条件================
SUBROUTINE INSCD(ND,NW,NFIX,LFIX,TK,F,DISP)
!ND总自由度,NW半带宽,NFIX总约束数,LFIX约束信息,TK总刚,F荷载列阵,DISP结点位移
	IMPLICIT NONE
	CHARACTER(LEN=50)::CTM="  4 IntroduceBoundaryConditions:"
    INTEGER::ND,NW,NFIX,I,J,II
    REAL*8,DIMENSION(:)::LFIX(NFIX),F(ND),DISP(ND)
    REAL*8,DIMENSION(:,:)::TK(ND,NW)
    DO I=1,NFIX
        II=LFIX(I)
        TK(II,1)=1.0
        F(II)=0.0
        DO J=II+1,II+NW-1
			TK(II,J-II+1)=0.0
        END DO
		DO J=II-1,II-NW+1,-1
			TK(J,II-J+1)=0.0
		END DO
    END DO
    DO I=1,ND
        DISP(I)=F(I)
    END DO
	PRINT '(A,10X,"DONE!")',CTM
	RETURN
END
!=================引入边界条件================

!===================求解方程==================
SUBROUTINE SOLVE(ND,NW,TK,F)	
	IMPLICIT NONE
	CHARACTER(LEN=50)::CTM="  5 Solve:"
	INTEGER::I,J,K,ND,NW
	REAL*8,DIMENSION(:,:)::TK(ND,NW)
	REAL*8,DIMENSION(:)::F(ND)
	DO K=1,ND-1
		DO I=K+1,MIN(K+NW-1,ND)
			F(I)=F(I)-((TK(K,I-K+1)*F(K))/TK(K,1))
			DO J=I,MIN(K+NW-1,ND)
				TK(I,J-I+1)=TK(I,J-I+1)-((TK(K,I-K+1)*TK(K,J-K+1))/TK(K,1))
				IF(ISNAN(TK(I,J-I+1)).OR.ISNAN(F(I))) THEN
					PRINT *,K,I,J-I+1,TK(K,I-K+1),F(K)
					GOTO 404
				ENDIF
			END DO
		END DO
	END DO
	F(ND)=F(ND)/TK(ND,1)
	DO I=ND-1,1,-1
		DO J=I+1,MIN(I+NW-1,ND)
			F(I)=F(I)-F(J)*TK(I,J-I+1)
		END DO
		F(I)=F(I)/TK(I,1)
		IF(ISNAN(F(I))) THEN
			PRINT *,I
			GOTO 404
		ENDIF
	END DO
	PRINT '(A,10X,"DONE!")',CTM
	RETURN
404	PRINT *,"ERROR"
	READ(*,*)
END
!=================求解方程===================

!=================单元应力===================
SUBROUTINE INNERFORCE(NN,NE,ND,E,ANU,T,GM,LOC,CX,CY,DISP)
!IE单元计数,NN结点总数,NE单元总数,ND总自由度,E弹性模量,ANU泊松比,T厚度,GM重度,LOC单元结点编码,CX,CY结点坐标,DISP结点位移
    IMPLICIT NONE
	CHARACTER(LEN=50)::CTM="  6 CalculateInerforce:"
    INTEGER::I,II   !I单元结点位移定位码,II总体位移定位码
    INTEGER::IE,NN,NE,ND
    REAL*8::E,ANU,T,GM,AE
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
        WRITE(80,301)IE,(INNF(I),I=1,3)
    END DO
301 FORMAT(3X,I4,4X,E13.3,2X,E13.3,2X,E13.3)
	PRINT '(A,10X,"DONE!")',CTM
	RETURN
END
!=================单元应力===================