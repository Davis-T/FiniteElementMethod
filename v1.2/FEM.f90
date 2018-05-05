MODULE CMD_PROGRESS
  IMPLICIT NONE
  PRIVATE
  LOGICAL , PARAMETER , PUBLIC :: CMD_PROGRESS_ABSOLUTE = .TRUE.
  TYPE , PUBLIC :: CLS_CMD_PROGRESS
    INTEGER , PRIVATE :: N , LENS , I
    CHARACTER :: M = "#" , O = "."
    CHARACTER(LEN=64) :: PREFIX
  CONTAINS
    PROCEDURE :: SET
    PROCEDURE :: PUT
END TYPE CLS_CMD_PROGRESS
  
CONTAINS

  SUBROUTINE SET( THIS , N , L )
    CLASS( CLS_CMD_PROGRESS ) :: THIS
    INTEGER , INTENT( IN ) :: N , L
    THIS % N    = N
    THIS % LENS = L
    THIS % I = 0
    THIS % PREFIX = " PROGRESS: " !//
  END SUBROUTINE SET
  
  SUBROUTINE PUT( THIS , K , BABSOL )
    CLASS( CLS_CMD_PROGRESS ) :: THIS
    INTEGER , INTENT( IN ) :: K
    LOGICAL , OPTIONAL :: BABSOL
    CHARACTER(LEN=1) :: BR
    INTEGER :: JM
    THIS % I = THIS % I + K
    IF ( PRESENT( BABSOL ) ) THEN
      IF ( BABSOL ) THIS % I = K
    END IF
    IF ( THIS % I > THIS % N ) THIS % I = THIS % N    
    JM = NINT( REAL( THIS%I * THIS%LENS ) / REAL( THIS%N ) )
    IF ( THIS%I < THIS%N ) THEN
      BR = CHAR(13)
    ELSE
      BR = CHAR(10)
    END IF
    WRITE( * , '(5A,F6.2,2A)',ADVANCE="NO") TRIM(THIS%PREFIX) , '[' , &
    !WRITE( * , '(5A,F6.2,2A\)') TRIM(THIS%PREFIX) , '[' , & !// 如您的编译器不支持，请用上方语句代替
      REPEAT(THIS%M , JM ) , REPEAT( THIS%O , THIS%LENS-JM ) , '] ' , THIS%I*100.0/THIS%N , "%" , BR
  END SUBROUTINE PUT
  
END MODULE CMD_PROGRESS
PROGRAM Finite_element_Analysis
    IMPLICIT NONE
    INTEGER::I
    INTEGER::NN,NE,ND,NFIX,NFORCE,NTYPE!NN结点总数,NE单元总数,ND总自由度数,NFIX约束数,NFORCE荷载总数,NTYPE求解类型（1平面应力问题、2平面应变问题）
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
    ALLOCATE(CX(NN),CY(NN),LFIX(NFIX),LOC(NE,3),TK(ND,ND),FN(NFORCE*2),F(ND),DISP(ND))
    CALL DATAIN(NN,NE,NFIX,NFORCE,LOC,CX,CY,LFIX,FN)   !读取结构信息——LOC单元结点编码、CX,CY结点坐标、LFIX约束信息、FN荷载信息
    CALL GLOSTIF(NN,NE,ND,E,ANU,T,LOC,CX,CY,TK)    !计算总刚矩阵TK
    CALL GLOFORC(NN,NE,ND,NFORCE,T,GM,LOC,CX,CY,FN,F)    !计算结点荷载列阵F
    CALL INSCD(ND,NFIX,LFIX,TK,F,DISP)     !引入边界条件
    CALL SOLVE(ND,TK,DISP)     !求解位移DISP
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
	WRITE(60,'(//,4X,"TotleTime:",F12.3,2X,"seconds")') TOTLETIME
	CLOSE(50)
    CLOSE(60)
	CLOSE(70)
	CLOSE(80)
	PRINT '(/,4X,A)', "AnalysisComplete!"
	READ(*,*)
END PROGRAM Finite_element_Analysis

!==================输入数据==================
SUBROUTINE DATAIN(NN,NE,NFIX,NFORCE,LOC,CX,CY,LFIX,FN)  
!NN结点总数,NE单元总数,NFIX约束数,NFORCE荷载总数,LOC单元编码,CX,CY编码坐标,LFIX约束信息,FN荷载信息
    IMPLICIT NONE
	CHARACTER(LEN=50)::CTM
    INTEGER::I,J,K,L
    INTEGER::NN,NE,NFIX,NFORCE
    REAL*8,DIMENSION(:)::CX(NN),CY(NN),LFIX(NFIX),FN(NFORCE)
    REAL*8,DIMENSION(:,:)::LOC(NE,3)
    READ(50,*)(LOC(I,1),LOC(I,2),LOC(I,3),I=1,NE)    !读取结点编码
    READ(50,*)(CX(J),CY(J),J=1,NN)       !读取结点坐标
    READ(50,*)(LFIX(K),K=1,NFIX)     !读取约束信息
    READ(50,*)(FN(L),L=1,NFORCE*2)       !读取荷载信息
	CTM="  1 ReadData:"
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
SUBROUTINE GLOSTIF(NN,NE,ND,E,ANU,T,LOC,CX,CY,TK)
!NN结点总数,NE单元总数,ND总自由度,E弹性模量,ANU泊松比,T厚度,LOC单元结点编码,CX,CY节点坐标,TK总刚
    IMPLICIT NONE
	CHARACTER(LEN=50)::CTM
    INTEGER::I,J,II,JJ    !II,JJ定位编码
    INTEGER::IE,NN,NE,ND    !IE单元计数,NN结点总数,NE单元总数,ND总自由度
    REAL*8::E,ANU,T   !E弹性模量,ANU泊松比,T单元结构厚度
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
                TK(II,JJ)=TK(II,JJ)+KE(I,J)
            END DO
        END DO
    END DO
	CTM="  2 MakeGlobalStiffnessMatrix:"
	PRINT '(A,10X,"DONE!")',CTM
	RETURN
END
!===============集成总体刚度矩阵===============

!=================计算结点荷载================
SUBROUTINE GLOFORC(NN,NE,ND,NFORCE,T,GM,LOC,CX,CY,FN,F)    
!NN结点总数,NE单元总数,ND自由度总数,NFORCE荷载总数,T厚度,GM重度,LOC单元结点编码,CX,CY结点坐标,FG等效重力荷载,FN荷载信息,F总结点荷载列阵
	IMPLICIT NONE
	CHARACTER(LEN=50)::CTM
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
	CTM="  3 MakeLoadArray:"
	PRINT '(A,10X,"DONE!")',CTM
	RETURN
END
!=================计算结点荷载================

!=================引入边界条件================
SUBROUTINE INSCD(ND,NFIX,LFIX,TK,F,DISP)
!ND总自由度,NFIX总约束数,LFIX约束信息,TK总刚,F荷载列阵,DISP结点位移
	IMPLICIT NONE
	CHARACTER(LEN=50)::CTM
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
    END DO
	CTM="  4 IntroduceBoundaryConditions:"
	PRINT '(A,10X,"DONE!")',CTM
	RETURN
END
!=================引入边界条件================

!==============Gauss_Jordan消去法求解方程==============
SUBROUTINE SOLVE(N,A,B)
	USE CMD_PROGRESS	
	IMPLICIT NONE
	CHARACTER(LEN=50)::CTM="    Gauss-Jordan elimination method:"
	INTEGER::I,J,K,N,IS
	INTEGER,DIMENSION(:)::JS(N)
	REAL*8,DIMENSION(:,:)::A(N,N)
	REAL*8,DIMENSION(:)::B(N)
	REAL*8::D,T
	TYPE( CLS_CMD_PROGRESS ) ::PROGRESS
	CALL PROGRESS % SET( N = N , L = 50 )
	PROGRESS % PREFIX = "  5 Solve:"
	DO K=1,N
		CALL PROGRESS % PUT( K , CMD_PROGRESS_ABSOLUTE )
		D=0.0
		DO I=K,N
			DO J=K,N
				IF (ABS(A(I,J)).GT.D) THEN
					D=ABS(A(I,J))
					JS(K)=J
					IS=I
				END IF
			END DO
		END DO
		DO J=K,N
			T=A(K,J)
			A(K,J)=A(IS,J)
			A(IS,J)=T
		END DO
		T=B(K)
		B(K)=B(IS)
		B(IS)=T
		DO I=1,N
			T=A(I,K)
			A(I,K)=A(I,JS(K))
			A(I,JS(K))=T
		END DO
		T=A(K,K)
		DO J=K+1,N
			IF (A(K,J).NE.0.0) THEN
				A(K,J)=A(K,J)/T
			ENDIF
		END DO
		B(K)=B(K)/T
		DO J=K+1,N
			IF (A(K,J).NE.0.0) THEN
				DO I=1,N
					IF ((I.NE.K).AND.(A(I,K).NE.0.0)) THEN
						  A(I,J)=A(I,J)-A(I,K)*A(K,J)
					END IF
				END DO
			END IF
		END DO
		DO I=1,N
			IF ((I.NE.K).AND.(A(I,K).NE.0.0)) THEN
				B(I)=B(I)-A(I,K)*B(K)
			END IF
		END DO
	END DO
	DO K=N,1,-1
		IF (K.NE.JS(K)) THEN
			T=B(K)
			B(K)=B(JS(K))
			B(JS(K))=T
		END IF
	END DO
	PRINT '(A,10X,"DONE!")',CTM
	RETURN
END
!==============Gauss_Jordan消去法求解方程==============

!=================单元应力===================
SUBROUTINE INNERFORCE(NN,NE,ND,E,ANU,T,GM,LOC,CX,CY,DISP)
!IE单元计数,NN结点总数,NE单元总数,ND总自由度,E弹性模量,ANU泊松比,T厚度,GM重度,LOC单元结点编码,CX,CY结点坐标,DISP结点位移
    IMPLICIT NONE
	CHARACTER(LEN=50)::CTM
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
	CTM="  6 CalculateInerforce:"
	PRINT '(A,10X,"DONE!")',CTM
	RETURN
END
!=================单元应力===================