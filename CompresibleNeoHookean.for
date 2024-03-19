C     Compresible NeoHookean UMAT
C     author = Sina Taghizadeh
C     license = GNU GPL Version 3.0
          
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
	
	PARAMETER(Zero=0.D0, One=1.D0, Two=2.D0, Three=3.D0, Nine=9.D0)
	
	DIMENSION DistortionGrad(NDI,NDI), BBar(NDI,NDI),
     2 StressTensor(NDI,NDI), Eye(NDI,NDI),
     1 JacobianMatrix(NDI,NDI,NDI,NDI), Index1(NTENS), Index2(NTENS)
	
	C10 = PROPS(1) ! First mechanical constant
	D1  = PROPS(2) ! Second mechanical constant
	
      ! Define identity matrix
	do i = 1, NDI
          do j = 1, NDI
              Eye(i,j) = Zero  ! Initialize all elements to zero
          end do
          Eye(i,i) = One  ! Set diagonal elements to one
      end do
      
      ! Define jacobin of deformation gradient
	 DetDefGrad = DFGRD1(1,1) * DFGRD1(2,2) * DFGRD1(3,3)
     1	        - DFGRD1(1,2) * DFGRD1(2,1) * DFGRD1(3,3)
     2            + DFGRD1(1,2) * DFGRD1(2,3) * DFGRD1(3,1)
     3	        + DFGRD1(1,3) * DFGRD1(3,2) * DFGRD1(2,1)
     4	        - DFGRD1(1,3) * DFGRD1(3,1) * DFGRD1(2,2)
     5	        - DFGRD1(2,3) * DFGRD1(3,2) * DFGRD1(1,1)
       
	! Define the distortion gradient
      Scale = (DetDefGrad)**(-One/Three)
      DO i = 1, NDI
        DO j = 1, NDI
          DistortionGrad(j,i) = Scale * DFGRD1(j,i)
        END DO
      END DO	
      
      ! Define the deviatoric left Cauchy-Green deformation tensor
	BBar = MATMUL(DistortionGrad, TRANSPOSE(DistortionGrad)) 
      
      ! Trace of The deviatoric left Cauchy-Green deformation tensor
	TraceBBar = BBar(1,1) + BBar(2,2) + BBar(3,3)  
	
      ! Define Stress tensor 3*3
	StressTensor = (Two*C10/DetDefGrad)*(BBar-(TraceBBar/Three)*Eye)
     1 +(Two/D1)*(DetDefGrad-One)*Eye 
      
      ! Define Abaqus stress vector 6*1
	STRESS(1) = StressTensor(1,1)
	STRESS(2) = StressTensor(2,2)
	STRESS(3) = StressTensor(3,3)
	STRESS(4) = StressTensor(1,2)
	STRESS(5) = StressTensor(1,3)
	STRESS(6) = StressTensor(2,3)
	
      ! Define system jacobian 3*3*3*3 matrix 
	DO i=1, NDI
		DO j=1, NDI
			DO k=1, NDI
				DO l=1, NDI
	JacobianMatrix(i,j,k,l) = (Two*C10/DetDefGrad)*((Eye(i,k)*BBar(j,l)
	1 +BBar(i,k)*Eye(j,l)+Eye(i,l)*BBar(j,k)+BBar(i,l)*Eye(j,k))/Two
	2 -(Two/Three)*Eye(i,j)*BBar(k,l)-(Two/Three)*BBar(i,j)*Eye(k,l)
	3 +(Two/Nine)*Eye(i,j)*Eye(k,l)*TraceBBar)
	4 +(Two/D1)*(Two*DetDefGrad-One)*Eye(i,j)*Eye(k,l)
				END DO
			END DO
		END DO
      END DO
	
      ! Define system jacobian 6*6 matrix or DDSDDE
	Index1 = [1, 2, 3, 1, 1, 2]
	Index2 = [1, 2, 3, 2, 3, 3]
	DO i=1, NTENS
		DO j=1, NTENS
		DDSDDE(i,j) = JacobianMatrix(Index1(i),Index2(i),Index1(j),Index2(j))
		END DO
	END DO
	
      RETURN
      END