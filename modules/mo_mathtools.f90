MODULE MO_MATHTOOLS
!! This module contains all the mathematical functions or subroutines for practical use
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE LAPACK95

IMPLICIT NONE

! Geo Constant
real*8, parameter :: pi=3.141592653589793d0
real*8, parameter :: RieFac = (1.0d0)/sqrt(2.0d0)  !!!! riemann transform factor

real*8, parameter :: AbsZeroTemp = -273.15d0 
real*8, parameter :: EarthG = 9.80d0
real*8, parameter :: EarthR = 6371.0d0 * 1000.0d0
real*8, parameter :: EarthOmega = 2.0d0*pi/(86400.0d0)
real*8, parameter :: EarthBeta  = 2.0d0*EarthOmega/EarthR
real*8, parameter :: OneDegLen = (EarthR)*(2.0d0*pi/360.0d0)

integer*8,parameter  :: Norm_DayMonSum (12) = (/0,31,59,90,120,151,181,212,243,273,304,334/)
integer*8,parameter  :: Leap_DayMonSum (12) = (/0,31,60,91,121,152,182,213,244,274,305,335/)
integer*8, parameter :: Norm_DayMon(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
integer*8, parameter :: Leap_DayMon(12) = (/31,29,31,30,31,30,31,31,30,31,30,31/)

!!! 31*2 character
character(len=62),parameter :: DateChar = &
                               '01'//'02'//'03'//'04'//'05'//'06'//'07'//'08'//'09'//'10'//&
                               '11'//'12'//'13'//'14'//'15'//'16'//'17'//'18'//'19'//'20'//&
                               '21'//'22'//'23'//'24'//'25'//'26'//'27'//'28'//'29'//'30'//'31'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! OVERLOAD
    INTERFACE Num2Str
        
        module procedure Int8_2Str
        module procedure Re8_2Str
        module procedure Int4_2Str
        module procedure Re4_2Str
        
    END INTERFACE Num2Str
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTERFACE Update_MonMe_Var
        
        module procedure Update_MonMe_Var0D
        module procedure Update_MonMe_Var1D
        module procedure Update_MonMe_Var2D
        module procedure Update_MonMe_Var3D
        
    END INTERFACE Update_MonMe_Var
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    CONTAINS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Update_MonMe_Var0D (Var, MonMe_Var, &
               Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Use_LeapYear, CYear, CMonth, Day_StepNum)
    
        real*8 :: Var,MonMe_Var
        
        integer*8 :: Arrive_CMonth_FLAG !!! Arrive_FLAG: arrive at the start of month
        integer*8 :: Reach_CMonth_FLAG  !!! Reach_FLAG:  arrive at the end of month
        
        integer*8 :: Use_LeapYear   !!! use leap year (1) or not (0)
        integer*8 :: CYear          !!! calendar year
        integer*8 :: CMonth         !!! calendar month
        
        integer*8 :: Day_StepNum    !!! Integration Step Num Each Day (e.g. Day_StepNum = Day_NT = 8)
        
        integer*8 :: Total_Step
        
        Total_Step = ( Day_StepNum * ( Norm_DayMon(CMonth) + Use_LeapYear * Math_IsLeapYear(CYear) ) );
        
        if ( Arrive_CMonth_FLAG == 1 ) then
            MonMe_Var = 0.0d0;
        end if
        
        MonMe_Var = MonMe_Var + Var;
        
        if ( Reach_CMonth_FLAG  == 1 ) then
             MonMe_Var = ( MonMe_Var )/dble( Total_Step );
        end if
    
    END SUBROUTINE Update_MonMe_Var0D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Update_MonMe_Var1D (Var, MonMe_Var, &
               Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Use_LeapYear, CYear, CMonth, Day_StepNum)
    
        real*8,allocatable :: Var(:),MonMe_Var(:)
        
        integer*8 :: Arrive_CMonth_FLAG !!! Arrive_FLAG: arrive at the start of month
        integer*8 :: Reach_CMonth_FLAG  !!! Reach_FLAG:  arrive at the end of month
        
        integer*8 :: Use_LeapYear   !!! use leap year (1) or not (0)
        integer*8 :: CYear          !!! calendar year
        integer*8 :: CMonth         !!! calendar month
        
        integer*8 :: Day_StepNum    !!! Integration Step Num Each Day (e.g. Day_StepNum = Day_NT = 8)
        
        integer*8 :: Total_Step
        
        Total_Step = ( Day_StepNum * ( Norm_DayMon(CMonth) + Use_LeapYear * Math_IsLeapYear(CYear) ) );
        
        if ( Arrive_CMonth_FLAG == 1 ) then
            MonMe_Var = 0.0d0;
        end if
        
        MonMe_Var = MonMe_Var + Var;
        
        if ( Reach_CMonth_FLAG  == 1 ) then
             MonMe_Var = ( MonMe_Var )/dble( Total_Step );
        end if
    
    END SUBROUTINE Update_MonMe_Var1D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Update_MonMe_Var2D (Var, MonMe_Var, &
               Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Use_LeapYear, CYear, CMonth, Day_StepNum)
    
        real*8,allocatable :: Var(:,:),MonMe_Var(:,:)
        
        integer*8 :: Arrive_CMonth_FLAG !!! Arrive_FLAG: arrive at the start of month
        integer*8 :: Reach_CMonth_FLAG  !!! Reach_FLAG:  arrive at the end of month
        
        integer*8 :: Use_LeapYear   !!! use leap year (1) or not (0)
        integer*8 :: CYear          !!! calendar year
        integer*8 :: CMonth         !!! calendar month
        
        integer*8 :: Day_StepNum    !!! Integration Step Num Each Day (e.g. Day_StepNum = Day_NT = 8)
        
        integer*8 :: Total_Step
        
        Total_Step = ( Day_StepNum * ( Norm_DayMon(CMonth) + Use_LeapYear * Math_IsLeapYear(CYear) ) );
        
        if ( Arrive_CMonth_FLAG == 1 ) then
            MonMe_Var = 0.0d0;
        end if
        
        MonMe_Var = MonMe_Var + Var;
        
        if ( Reach_CMonth_FLAG  == 1 ) then
             MonMe_Var = ( MonMe_Var )/dble( Total_Step );
        end if
    
    END SUBROUTINE Update_MonMe_Var2D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Update_MonMe_Var3D (Var, MonMe_Var, &
               Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Use_LeapYear, CYear, CMonth, Day_StepNum)
    
        real*8,allocatable :: Var(:,:,:),MonMe_Var(:,:,:)
        
        integer*8 :: Arrive_CMonth_FLAG !!! Arrive_FLAG: arrive at the start of month
        integer*8 :: Reach_CMonth_FLAG  !!! Reach_FLAG:  arrive at the end of month
        
        integer*8 :: Use_LeapYear   !!! use leap year (1) or not (0)
        integer*8 :: CYear          !!! calendar year
        integer*8 :: CMonth         !!! calendar month
        
        integer*8 :: Day_StepNum    !!! Integration Step Num Each Day (e.g. Day_StepNum = Day_NT = 8)
        
        integer*8 :: Total_Step
        
        Total_Step = ( Day_StepNum * ( Norm_DayMon(CMonth) + Use_LeapYear * Math_IsLeapYear(CYear) ) );
        
        if ( Arrive_CMonth_FLAG == 1 ) then
            MonMe_Var = 0.0d0;
        end if
        
        MonMe_Var = MonMe_Var + Var;
        
        if ( Reach_CMonth_FLAG  == 1 ) then
             MonMe_Var = ( MonMe_Var )/dble( Total_Step );
        end if
    
    END SUBROUTINE Update_MonMe_Var3D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function Int8_2Str(Num) result(Str)
    
        integer*8 :: Num
        character(len=100) :: Str
        
        write(Str,'(I8)'),Num
        Str = AdjustL(trim(Str))
    
    end function Int8_2Str
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function Re8_2Str(Num) result(Str)
    
        real*8 :: Num
        character(len=100) :: Str
        
        write(Str,'(f13.5)'),Num
        Str = AdjustL(trim(Str))
        
        IF(Num >= 0.0d0 )THEN
        
           Str = ' '//Str;
        
        END IF
    
    end function Re8_2Str
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function Int4_2Str(Num) result(Str)
    
        integer*4 :: Num
        character(len=100) :: Str
        
        write(Str,'(I8)'),Num
        Str = AdjustL(trim(Str))
    
    end function Int4_2Str
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function Re4_2Str(Num) result(Str)
    
        real*4 :: Num
        character(len=100) :: Str
        
        write(Str,'(f13.5)'),Num
        Str = AdjustL(trim(Str))
        
        IF(Num >= REAL(0.0d0) )THEN
        
           Str = ' '//Str;
        
        END IF
    
    end function Re4_2Str
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Math_YMD_To_Date(CYear,CMonth,CDay) RESULT(DateStr)
    
        !!! calendar Year,Month,Day
        integer*8 :: CYear
        integer*8 :: CMonth
        integer*8 :: CDay
        
        character(len=100) :: DateStr
        character(len=4)   :: YearStr
        
        write(YearStr, '(I4.4)') CYear
        
        DateStr = &
        trim(YearStr)//'-'//DateChar(2*CMonth-1:2*CMonth)//'-'//DateChar(2*CDay-1:2*CDay)
        
        DateStr = (trim(DateStr));
    
    END FUNCTION Math_YMD_To_Date
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Math_IsLeapYear(y) RESULT (FLAG)
    integer*8 :: y
    integer*8 :: FLAG
    
        if (mod(y, 4) == 0) then
            if (mod(y, 100) /= 0 .or. mod(y, 400) == 0) then
                FLAG = 1
            else
                FLAG = 0
            end if
        else
            FLAG = 0
        end if
        
    END FUNCTION Math_IsLeapYear
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Math_TryFunc(NY,NX,FuncName) RESULT(TFunc)
    !!!! generate widely used trial function
         integer*8 :: NY,NX  !!! TFunc has (NX+1) elements
         character :: FuncName !!! 'G': 2D gauss distribution
         real*8,allocatable :: TXFunc(:,:),TYFunc(:,:),TFunc(:,:) !!! (1,NX+1), be careful
         real*8,allocatable :: Xspac(:,:),Yspac(:,:) !!! (1,NX+1), be careful
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         allocate(Xspac(1,NX+1))
         allocate(Yspac(NY+1,1))
         
         allocate(TXFunc(1,NX+1),TYFunc(NY+1,1))
         allocate(TFunc(NY+1,NX+1))
         
         Xspac(1,1:NX+1) = Math_Linspace(0.0d0, 1.0d0, NX+1)
         Yspac(1:NY+1,1) = Math_Linspace(1.0d0,-1.0d0, NY+1)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if(FuncName == 'G') then
         
           TXFunc = exp( - (20.0d0) * (Xspac - ( (0.5d0)* Xspac(1,1)+ (0.5d0)* Xspac(1,NX+1) ) )**2 )
           TYFunc = exp( - (20.0d0) * (Yspac - ( (0.5d0)* Yspac(1,1)+ (0.5d0)* Yspac(NY+1,1) ) )**2)
           
           TFunc = matmul(TYFunc,TXFunc)
         
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!! DEALLOCATE
         deallocate(Xspac,Yspac)
         deallocate(TXFunc,TYFunc)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    END FUNCTION Math_TryFunc
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! index check
    !!! it indicates natural(non-gradient) boundary condition outside
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION SAFE(IND,NIND) RESULT(SAFEIND)
    
         integer*8 :: IND,NIND,SAFEIND
         
         SAFEIND = max(IND,1);
         SAFEIND = min(SAFEIND,NIND);
    
    END FUNCTION SAFE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! smabs
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION smabs(X) RESULT(Y)
    !!! smooth absolute function to make abs(x) differentiable 
    
         real*8 :: X,Y
         Y = sqrt( X*X + (1.0d-15)*(1.0d-15) );
    
    END FUNCTION smabs
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! smsign
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION smsign(X) RESULT(Y)
    !!! smooth sign function to make sign(x) differentiable
    
        real*8 :: X,Y
        Y = 2.0d0/(1.0d0 + exp(-1.0d2*X))-1.0d0;
    
    END FUNCTION smsign
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! smheav
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION smheav(X) RESULT(Y)
    !!! smooth heavside function to make H(x) differentiable
    
        real*8 :: X,Y
        Y = 1.0d0/(1.0d0 + exp(-1.0d2*X))
    
    END FUNCTION smheav
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! smmax
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION smmax(X,Y) RESULT(Z)
    !!! smooth max(X,Y)
    
        real*8 :: X,Y,Z
        Z = smheav(X-Y)*X + smheav(Y-X)*Y
    
    END FUNCTION smmax
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! smmin
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION smmin(X,Y) RESULT(Z)
    !!! smooth min(X,Y)
    
        real*8 :: X,Y,Z
        Z = smheav(X-Y)*Y + smheav(Y-X)*X
        
    END FUNCTION smmin
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! psgn_skew
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    FUNCTION psgn_skew (y,a,b,c,m) RESULT(z)

        real*8 :: y,a,b,c,m
        real*8 :: z

        z =  b*(+tanh(a/b) - tanh((a-y)/b));
        z = z - m;
        z = c * z;

    END FUNCTION psgn_skew
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! msgn_skew 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION msgn_skew (y,a,b,c,m) RESULT(z)

        real*8 :: y,a,b,c,m
        real*8 :: z

        z =  b*(-tanh(a/b) + tanh((a+y)/b));
        z = z - m;
        z = c * z;

    END FUNCTION  msgn_skew
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! psgn
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION psgn(X) RESULT(Y)
        
          real*8 :: X,Y
          !!! Y = 0.5d0*( X + abs(X) );
          Y = 0.5d0*( X + abs(X) );
    
    END FUNCTION psgn
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! msgn
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION msgn(X) RESULT(Y)
          
          real*8 :: X,Y
          !!! Y = 0.5d0*( X - abs(X) );
          Y = 0.5d0*( X - abs(X) );
    
    END FUNCTION msgn
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! UbdCtrl: upper bound control
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION UbdCtrl(X,Ubd) RESULT(Y)
    
        !!! UbdCtrl(X,0.0d0) = msgn(X)
        
        real*8 :: X,Y,Ubd
        !!! Y = Ubd - (0.5d0)*( (Ubd-X) + abs(Ubd-X) );
        Y = Ubd - (0.5d0)*( (Ubd-X) + smabs(Ubd-X) );
    
    END FUNCTION UbdCtrl
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! LbdCtrl: lower bound control
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION LbdCtrl(X,Lbd) RESULT(Y)
    
        !!! LbdCtrl(X,0.0d0) = psgn(X)
        
        real*8 :: X,Y,Lbd
        !!! Y = Lbd + (0.5d0)*( (X-Lbd) + abs(X-Lbd) );
        Y = Lbd + (0.5d0)*( (X-Lbd) + smabs(X-Lbd) );
    
    END FUNCTION LbdCtrl
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Math_Eye
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Math_Eye(N) RESULT(EyeMat)
         implicit none
         integer*8 :: N,ind
         real*8,allocatable :: EyeMat(:,:)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         allocate(EyeMat(N,N));
         EyeMat = 0.0d0;
         do ind = 1,N
            EyeMat(ind,ind) = 1.0d0;
         end do
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END FUNCTION
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Math_Linspace
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Math_Linspace(LeftPoint,RightPoint,N) RESULT(Vec)
         implicit none
         real*8 :: LeftPoint  ! left end point
         real*8 :: RightPoint ! right end point
         integer*8 :: N       ! num of points
         real*8,allocatable :: Vec(:)  ! size : (N)
         integer*8 :: ind
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         allocate(Vec(N)) ! column major form
         
         if(N>=2)then
           Vec(1:N)=((RightPoint - LeftPoint)/(dble(N-1)))*dble((/(ind,ind=0,N-1)/)) + LeftPoint
         else
           Vec(1) = LeftPoint;
         end if
    
    END FUNCTION Math_Linspace
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Math_Powspace
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Math_Powspace(LeftPoint,RightPoint,N,Pow) RESULT(Vec)
    !!! to get the power distribution around the midpoint
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        real*8 :: LeftPoint, RightPoint, MidPoint
        integer*8 :: N
        real*8 :: Pow 
        real*8,allocatable :: Vec(:)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8,allocatable :: LeftHalfVec(:),RightHalfVec(:)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate(Vec(N), LeftHalfVec(N/2+1), RightHalfVec(N/2+1) )
        MidPoint = (0.5d0)*(LeftPoint + RightPoint)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        if(mod(N,2)==0)then
        
           RightHalfVec = MidPoint + &
           ( (RightPoint - MidPoint) )*( (Math_Linspace(0.0d0,1.0d0,N/2+1))**(dble(Pow)) );
           
           LeftHalfVec = MidPoint - & 
           ( (MidPoint - LeftPoint) )*( (Math_Linspace(0.0d0,1.0d0,N/2+1))**(dble(Pow)) );
           
           LeftHalfVec = flipVec(LeftHalfVec)
           
           Vec(1:N/2) = LeftHalfVec(1:N/2)
           
           Vec((N/2+1):N) = RightHalfVec(2:(N/2+1))
        
        else 
        
           RightHalfVec = MidPoint + &
           ( (RightPoint - MidPoint) )*( (Math_Linspace(0.0d0,1.0d0,N/2+1))**(dble(Pow)) );
           
           LeftHalfVec = MidPoint - & 
           ( (MidPoint - LeftPoint) )*( (Math_Linspace(0.0d0,1.0d0,N/2+1))**(dble(Pow)) );
           
           LeftHalfVec = flipVec(LeftHalfVec)
           
           Vec(1:N/2) = LeftHalfVec(1:N/2)
           
           Vec((N/2+1):N) = RightHalfVec(1:(N/2+1))
        
        end if
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        deallocate(LeftHalfVec,RightHalfVec)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    END FUNCTION Math_Powspace
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Math_Diff
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Math_Diff(vec) RESULT(diffvec)
    !!! just like diff in matlab
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        real*8,allocatable :: vec(:) 
        real*8,allocatable :: diffvec(:)
        integer*8 :: Dim
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Dim=size(vec)-1
        allocate(diffvec(Dim))
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        diffvec(1:Dim)=vec(2:Dim+1)-vec(1:Dim)
    END FUNCTION Math_Diff
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Math_InProd
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Math_InProd(FY,GY,Y) RESULT(InProd)
    !!! calculate the InProd of F,G
    !!! InProd = <F,G>/<G,G> or <F,G>
    !!! using the simplest trapezoidal integration method for equidistant nodes
    !!!
    !!! for Powspace distribution (see Math_Powspace)
    !!! we can use transform: Y = MidY + (HalfY)* Z^(P)
    !!! int{ F(Y)*G(Y)dY } = S{ F(Y)*G(Y)*[P*Zm^(P-1)]dZ }
    !!! however, its accuracy does not satisfy our needs
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        integer*8 :: NY,N
        real*8 :: MidY, HalfY, Pow = 2.0d0
        real*8,allocatable :: FY(:)  ! F(Y) (NY+1)
        real*8,allocatable :: GY(:)  ! G(Y) (NY+1)
        real*8,allocatable :: Y(:)   ! Y Nodes, (NY+1) Nodes
        
        real*8,allocatable :: Z(:)   ! Z Nodes, (NY+1) Nodes
        real*8,allocatable :: Zm(:)  ! mid-points of Z Nodes (NY+1)
        
        real*8,allocatable :: DiffY(:),DiffZ(:)
        real*8 :: FxG, GxG, FxGxZ,GxGxZ, InProd
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        NY=size(Y,dim=1)-1
        N = NY+1
        MidY = (0.5d0)*(Y(1) + Y(NY+1))
        HalfY = (0.5d0)*(Y(NY+1) - MidY)
        allocate(Z(NY+1),Zm(NY),DiffY(NY),DiffZ(NY))
        
        !!! transform Y to Z
        !if(mod(N,2)==0)then
        !
        !   Z(1:N/2) = (-1.0d0)*( (abs((Y(1:N/2)-MidY)/HalfY))**(1.0d0/Pow) )
        !   Z(N/2+1:N) = (+1.0d0)*( (abs((Y(N/2+1:N)-MidY)/HalfY))**(1.0d0/Pow) )
        !
        !else
        !
        !   Z(1:N/2) = (-1.0d0)*( (abs((Y(1:N/2)-MidY)/HalfY))**(1.0d0/Pow) )
        !   Z(N/2+1) = 0.0d0;
        !   Z(N/2+2:N) = (+1.0d0)*( (abs((Y(N/2+2:N)-MidY)/HalfY))**(1.0d0/Pow) )
        !
        !end if
        !
        !
        !Zm = (0.5d0)*(Z(1:N-1) + Z(2:N))
        !DiffZ  = abs(Math_Diff(Z))
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! simplest trapezoidal integration
        DiffY(1:NY)=abs(Math_Diff(Y))
        
        FxG = (0.5d0)*sum(DiffY * FY(1:NY)   * GY(1:NY)) + &
              (0.5d0)*sum(DiffY * FY(2:NY+1) * GY(2:NY+1))
        !FxG = (0.5d0)*FxG
        
        GxG = (0.5d0)*sum(DiffY * GY(1:NY)   * GY(1:NY)) + &
              (0.5d0)*sum(DiffY * GY(2:NY+1) * GY(2:NY+1))
        !GxG = (0.5d0)*GxG
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !FxGxZ = sum( (HalfY*Pow) * (abs(Zm)**(Pow-1.0d0)) *DiffZ * FY(1:NY) * GY(1:NY) ) + &
        !        sum( (HalfY*Pow) * (abs(Zm)**(Pow-1.0d0)) *DiffZ * FY(2:NY+1) * GY(2:NY+1) )
        !
        !FxGxZ = 0.5d0*FxGxZ
        !
        !GxGxZ = sum( (HalfY*Pow) * (abs(Zm)**(Pow-1.0d0)) *DiffZ * GY(1:NY) * GY(1:NY) ) + &
        !        sum( (HalfY*Pow) * (abs(Zm)**(Pow-1.0d0)) *DiffZ * GY(2:NY+1) * GY(2:NY+1) )
        !
        !GxGxZ = 0.5d0*GxGxZ
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        InProd = FxG
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        deallocate(Z,Zm,DiffY,DiffZ)
    
    END FUNCTION Math_InProd
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SPECIAL FUNCTIONS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Math_GenPCFunc
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Math_GenPCFunc(Y,M) RESULT(PCF_valmat)
    !!! generates the 0 to m order std parabolic cylinder function
    !!! the std PC function (PCF) follows the recursive relation
    !!! PCF(m,y)=sqrt(1/m)*sqrt(2)*y*PCF(m-1,y)-sqrt(m-1/m)*PCF(m-2,y)
    !!!
    !!! Using Schmidt orthogonalization to ensure orthogonality between Discrete Vectors
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    real*8,allocatable :: PCF_valmat(:,:)           ! PCF tmp val mat : (NY+1,M+1)
    integer*8  :: M                                 ! PCF order
    real*8,allocatable :: Y(:,:)                    ! NdY: (NY+1,1), descend order
    real*8,allocatable :: YVec(:)                   
    integer*8  :: NY,Mind,Nind,tmpMind     
    real*8,allocatable :: tmpAVec(:),tmpBVec(:),tmpCVec(:)   
    
    ! (NY+1),for Schmidt orthogonalization
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    NY=size(Y,dim=1)-1
    allocate(PCF_valmat(NY+1,M+1));
    allocate(YVec(NY+1),tmpAVec(NY+1),tmpBVec(NY+1),tmpCVec(NY+1));
    
    YVec = Y(1:NY+1,1)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(M>=2)then
        
      PCF_valmat(1:NY+1,1:1) = ((pi)**(-0.25d0))*exp(-(0.50d0)*Y*Y);
      PCF_valmat(1:NY+1,2:2) = ((pi)**(-0.25d0))*sqrt(2.0d0)*(Y*exp(-(0.50d0)*Y*Y));
      
      do Mind=2,M
          
         tmpMind=Mind+1;
         
         PCF_valmat(1:NY+1,tmpMind:tmpMind)= &
         
         sqrt(1/dble(Mind))*sqrt(2.0d0)*(Y*PCF_valmat(1:NY+1,tmpMind-1:tmpMind-1)) + &
         
         (-1)*sqrt(dble(Mind-1)/dble(Mind))*PCF_valmat(1:NY+1,tmpMind-2:tmpMind-2);
         
         
      end do
      
    elseif(M==1)then
        
        PCF_valmat(1:NY+1,1:1) = ((pi)**(-0.25d0))*exp(-(0.50d0)*Y*Y);
        PCF_valmat(1:NY+1,2:2) = ((pi)**(-0.25d0))*sqrt(2.0d0)*(Y*exp(-(0.50d0)*Y*Y));
        
    else ! M ==0
        
        PCF_valmat(:,1:1) = ((pi)**(-0.25d0))*exp(-(0.50d0)*Y*Y);
        
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Schmidt orthogonalization
    tmpAVec = PCF_valmat(1:NY+1,1)
    tmpBVec = tmpAVec
    
    !!! step one : orthogonalization
    if(M>=1)then
      do Mind = 1,M
         
         tmpAVec = PCF_valmat(1:NY+1,Mind+(1))
         tmpBVec = tmpAVec
         
         do Nind = 0,Mind-1
            
            tmpCVec = PCF_valmat(1:NY+1,Nind+(1))
            tmpBVec = tmpBVec - &
            (Math_InProd(tmpAVec,tmpCVec,YVec)/(Math_InProd(tmpCVec,tmpCVec,YVec)))*tmpCVec
         
         end do !! Mind
         
         PCF_valmat(1:NY+1,Mind+(1)) = tmpBVec
         
      end do
      
    end if !! Nind
    
    !!! step two : standardlization
    do Mind = 0,M
    
       tmpAVec = PCF_valmat(1:NY+1,Mind+(1))
       
       PCF_valmat(1:NY+1,Mind+(1)) = tmpAVec/sqrt(Math_InProd(tmpAVec,tmpAVec,YVec))
       
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! check orthogonalization
    !tmpAVec = PCF_valmat(1:NY+1,1)
    !tmpBVec = PCF_valmat(1:NY+1,1)
    !
    !write(*,*) Math_InProd(tmpAVec,tmpBVec,YVec)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    deallocate(YVec,tmpAVec,tmpBVec,tmpCVec)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    END FUNCTION Math_GenPCFunc
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! flipMatUD (UP and Down)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION flipMatUD(Mat) RESULT(FlipMat)
        implicit none
        real*8,allocatable :: Mat(:,:)
        real*8,allocatable :: FlipMat(:,:)
        integer*8 :: Dim01, Dim02, RIND
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Dim01 = size(Mat,dim=1);
        Dim02 = size(Mat,dim=2);
        allocate(FlipMat(Dim01,Dim02));
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO RIND = 1,Dim01
           FlipMat(RIND,1:Dim02)=Mat(Dim01-(RIND-1),1:Dim02);
        END DO
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END FUNCTION flipMatUD
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! flipMatLR (Left and Right)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION flipMatLR(Mat) RESULT(FlipMat)
        implicit none
        real*8,allocatable :: Mat(:,:)
        real*8,allocatable :: FlipMat(:,:)
        integer*8 :: Dim01, Dim02, CIND
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Dim01 = size(Mat,dim=1);
        Dim02 = size(Mat,dim=2);
        allocate(FlipMat(Dim01,Dim02));
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO CIND = 1,Dim02
           FlipMat(1:Dim01,CIND)=Mat(1:Dim01,Dim02-(CIND-1));
        END DO
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END FUNCTION flipMatLR
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! flipVec
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION flipVec(Vec) RESULT(FVec)
        implicit none
        real*8,allocatable :: Vec(:)
        real*8,allocatable :: FVec(:)
        integer*8          :: Dim, IND
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Dim = size(Vec)
        allocate(FVec(Dim))
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO IND = 1,Dim
           FVec(IND) = Vec(Dim-(IND-1));
        END DO
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END FUNCTION flipVec
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Bubble Sort 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Math_BubbleSort(arr,ind_arr,job) 
    !!!! return the right index
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         implicit none
         real*8,allocatable :: arr(:),arr_help(:)
         real*8 :: tmpval
         character :: job !!!! 'A' for "ascend", 'D' for "descend" , 'T' for "twin"
         integer*8,allocatable :: ind_arr(:),ind_arr_help(:)
         integer*8 :: N,i,j,tmpind
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         N = size(arr);
         allocate(arr_help(N)); arr_help = 0.0d0;
         allocate(ind_arr_help(N)); ind_arr_help = 0;
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ind_arr(1:N) = (/(i,i=1,N)/)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF(job == 'A')then
         
           DO i = 1,N
           
              DO j = 2, (N-i+1)
              
                 if(arr(j-1) > arr(j))then
                 
                   tmpval = arr(j-1);
                   arr(j-1) = arr(j);
                   arr(j) = tmpval;
                   
                   tmpind = ind_arr(j-1);
                   ind_arr(j-1) = ind_arr(j);
                   ind_arr(j) = tmpind;
                 
                 end if
                 
              END DO
              
           END DO
         
         END IF
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF(job == 'D')then
         
           DO i = 1,N
           
              DO j = 2, (N-i+1)
              
                 if(arr(j-1) <= arr(j))then
                 
                   tmpval = arr(j-1);
                   arr(j-1) = arr(j);
                   arr(j) = tmpval;
                   
                   tmpind = ind_arr(j-1);
                   ind_arr(j-1) = ind_arr(j);
                   ind_arr(j) = tmpind
                 
                 end if
                 
              END DO
              
           END DO
         
         END IF
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF(job == 'T')then
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!! step one : bubble sort with descend job
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DO i = 1,N
           
                DO j = 2, (N-i+1)
              
                if(arr(j-1) <= arr(j))then
                 
                    tmpval = arr(j-1);
                    arr(j-1) = arr(j);
                    arr(j) = tmpval;
                   
                    tmpind = ind_arr(j-1);
                    ind_arr(j-1) = ind_arr(j);
                    ind_arr(j) = tmpind
                 
                end if
                 
                END DO
              
            END DO
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!! step two : make twin value together in descending way
         !!! e.g. (+100,-100, + 50, -50, +20, -20 ,10, -10)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             DO i = 1,N/2
                arr_help(2*i-1) = arr(i);
                arr_help(2*i)   = arr(N-i+1);
            
                ind_arr_help(2*i-1) = ind_arr(i);
                ind_arr_help(2*i)   = ind_arr(N-i+1);
             END DO
         
             arr = arr_help;
             ind_arr = ind_arr_help
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         END IF 
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         deallocate(arr_help,ind_arr_help)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
    END SUBROUTINE Math_BubbleSort
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! FindCoord_In_Mseq
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Math_FindCoord_In_MFSeq(Zseq,MFseq,Value,Tone) RESULT(Coord)
    !!! find coord in monotone value sequence
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        real*8,allocatable :: Zseq(:) !!! coord 
        real*8,allocatable :: MFseq(:),Mseq(:) !!! monotone value sequence : MFseq
        !!! Mseq(i) <= Mseq(i+1)
        real*8 :: Val,Value,Coord,tmpK,tmpB
        integer*8 :: N,ind,flag
        character :: Tone             !!! 'A':ascending ,'D':descending for Fseq
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        N = size(Zseq)
        allocate(Mseq(N)); 
        
        if(Tone == 'A')then
           Mseq = MFseq; 
           Val = Value;
        elseif(Tone == 'D')then
           Mseq = (-1.0d0)*MFseq; 
           Val  = (-1.0d0)*Value;
        end if
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! tmpK = ( Mseq(i+1) - Mseq(i) )/( Zseq(i+1) - Zseq(i) )
        !!! tmpB = Mseq(i)
        !!! Y = tmpK*(Z - Zseq(i)) + tmpB
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!! linear extrapolation 
        if (Val < Mseq(1))then
           tmpK = ( Mseq(2) - Mseq(1) )/( Zseq(2) - Zseq(1) )
           tmpB = Mseq(1)
           
           Coord = (Val - tmpB)/(tmpK) + Zseq(1)
           
           return
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (Val >= Mseq(N))then
           tmpK = ( Mseq(N) - Mseq(N-1))/( Zseq(N) - Zseq(N-1) )
           tmpB = Mseq(N-1)
           
           Coord = (Val - tmpB)/(tmpK) + Zseq(N-1)
           return
           
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!! linear interpolation
        flag = 0;
        ind = 1; 
        
        DO WHILE (flag .EQ. 0)
           if( (Val >= Mseq(ind)) .AND. (Val < Mseq(ind+1)) )then
               
               tmpK = ( Mseq(ind+1) - Mseq(ind) )/( Zseq(ind+1) - Zseq(ind) )
               tmpB = Mseq(ind)
           
               Coord = (Val - tmpB)/(tmpK)  + Zseq(ind)
               
               flag = 1;
               
           end if
           
           ind = ind +1;
           
        END DO
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
    END FUNCTION Math_FindCoord_In_MFseq
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! FindValue_In_MZSeq
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Math_FindValue_In_MZSeq(MZseq,Fseq,Coord,Tone) RESULT(Value)
    !!! find value in montone coord sequence (MZseq) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        real*8,allocatable :: Zseq(:),MZseq(:) !!! montone coord sqeuence : MZseq
        !!! Zseq(i) <= Zseq(i+1)
        real*8,allocatable :: Fseq(:)          !!! value sequence
        real*8 :: Value,Coord,Coo,tmpK,tmpB
        integer*8 :: N,ind,flag
        character :: Tone             !!! 'A':ascending ,'D':descending for Zseq
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        N = size(Fseq)
        allocate(Zseq(N));
        
        if(Tone == 'A')then
        
           Zseq = MZseq;
           Coo  = Coord
           
        elseif(Tone == 'D')then
        
           Zseq = (-1.0d0)*MZseq;
           Coo  = (-1.0d0)*Coord;
        
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! tmpK = ( Fseq(i+1) - Fseq(i) )/( Zseq(i+1) - Zseq(i) )
        !!! tmpB = Fseq(i)
        !!! Y = tmpK*(Z - Zseq(i)) + tmpB
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!! linear extrapolation 
        if (Coo < Zseq(1))then
        
           tmpK = ( Fseq(2) - Fseq(1) )/( Zseq(2) - Zseq(1) )
           tmpB = Fseq(1)
           
           Value = tmpK *(Coo - Zseq(1)) + tmpB
           
           !!! be careful (natural boundary condition)
           Value = Fseq(1)
           
           !!! be careful (zero boundary condition)
           !!! Value = 0.0d0;
           
           return
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elseif (Coo >= Zseq(N))then
           
           tmpK = ( Fseq(N) - Fseq(N-1))/( Zseq(N) - Zseq(N-1) )
           tmpB = Fseq(N-1)
           
           Value = tmpK *(Coo - Zseq(N-1)) + tmpB
           
           !!! be careful (natural boundary condition)
           Value = Fseq(N)
           
           !!! be careful (zero boundary condition)
           !!! Value = 0.0d0;
           
           return
           
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else
        
            !!!!!! linear interpolation
            flag = 0;
            ind = 1; 
        
            DO WHILE ( (flag .EQ. 0).AND.(ind<=(N-1)) )
               if( (Coo >= Zseq(ind)) .AND. (Coo < Zseq(ind+1)) )then
               
                   tmpK = ( Fseq(ind+1) - Fseq(ind))/( Zseq(ind+1) - Zseq(ind) )
                   tmpB = Fseq(ind)
           
                   Value = tmpK *(Coo - Zseq(ind)) + tmpB
               
               
                   flag = 1;
               
               end if
           
               ind = ind + 1;
           
             END DO
             
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    END FUNCTION Math_FindValue_In_MZSeq
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Math_RegTriDiagSolver
    FUNCTION Math_RegTriDiagSolver(AIN,BIN) RESULT(X)
    !!! colomn principal element elimination method to solver regular tri-diag matrix
        implicit none
        integer*8 :: j,i,N  ! N must be >= 3
        real*8    :: swaptmp,tmpQ
        real*8,allocatable :: A(:,:),AIN(:,:),B(:),BIN(:),X(:) 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        N = size(AIN,dim=1)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (N<=2)then
           write(*,*) " N must be >= 3, Math_RegTriDiagSolver stopped "
           return
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        allocate(X(N));   X=0.0d0
        allocate(A(N,N)); A=AIN;
        allocate(B(N));   B=BIN;
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO i = 1, N-1
           !!! select colomn principal element   
           if(abs(A(i+1,i))>abs(A(i,i)))then
      
            swaptmp = A(i,i); A(i,i)=A(i+1,i); A(i+1,i)=swaptmp;
        
            swaptmp = A(i,i+1); A(i,i+1)=A(i+1,i+1); A(i+1,i+1)=swaptmp;
        
            if(i<=(N-2))then
            
              swaptmp = A(i,i+2); A(i,i+2)=A(i+1,i+2); A(i+1,i+2) = swaptmp; 
              
            end if
            
            swaptmp = B(i); B(i)=B(i+1); B(i+1)=swaptmp;
        
            end if
            
            !!! elimination
            tmpQ = A(i+1,i)/A(i,i);
    
            A(i+1,i)   =   A(i+1,i) - A(i,i)*tmpQ;
            A(i+1,i+1) =   A(i+1,i+1) - A(i,i+1)*tmpQ;
    
            if(i<=(N-2))then 
               A(i+1,i+2) = A(i+1,i+2) - A(i,i+2)*tmpQ;
            end if
    
            B(i+1) = B(i+1) - B(i)*tmpQ;
    
        END DO
       
        !!! go back
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        X(N)   = B(N)/A(N,N);
        X(N-1) = (B(N-1)-X(N)*A(N-1,N))/A(N-1,N-1);

        DO i=(N-2),1,(-1)
        
           X(i)=(B(i) - X(i+1)*A(i,i+1) - X(i+2)*A(i,i+2))/A(i,i);

        END DO
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        deallocate(A,B)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    END FUNCTION Math_RegTriDiagSolver
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Math_CycTriDiagSolver
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Math_CycTriDiagSolver(A,B) RESULT(X)
    !!! colomn principal element elimination method to solver cyclic tri-diag matrix
        implicit none
        real*8,allocatable :: X(:),Z(:),Y(:),SUBA(:,:),B1(:),B2(:)
        real*8,allocatable :: A(:,:),B(:)
        integer*8 :: N
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        N = size(A,dim=1);
        allocate(SUBA(N-1,N-1),X(N),Z(N-1),Y(N-1),B1(N-1),B2(N-1));
        SUBA=0.0d0; X=0.0d0; Z=0.0d0; Y=0.0d0; B1=0.0d0; B2=0.0d0;
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBA(1:N-1,1:N-1)=A(1:N-1,1:N-1);
        
        B1(1)=(-A(1,N)); 
        B1(N-1)=(-A(N-1,N));
        B2(1:N-1)=B(1:N-1);
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Z = Math_RegTriDiagSolver(SUBA,B1);

        Y = Math_RegTriDiagSolver(SUBA,B2);
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        X(N) =(B(N)-(A(N,1)*Y(1) + A(N,N-1)*Y(N-1)))/(A(N,N) +(A(N,1)*Z(1) + A(N,N-1)*Z(N-1))); 

        X(1:N-1) = X(N)*Z(1:N-1) + Y(1:N-1);
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        deallocate(SUBA,Z,Y,B1,B2)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    END FUNCTION Math_CycTriDiagSolver
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! LAPACK ZGEEV
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Math_Eig(A,RVEC,LAM)
        !!! Compute the eigvenvalue and eigenvector of double precision matrix A
        !!! Using ZGEEV
        !!! See http://www.math.utah.edu/software/lapack/lapack-z/zgeev.html
        
        !NAME
        !   ZGEEV - compute for an N-by-N complex nonsymmetric matrix A,
        !   the eigenvalues and, optionally, the left and/or right
        !   eigenvectors
        !
        ! SYNOPSIS
        !      SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR,
        !                        LDVR, WORK, LWORK, RWORK, INFO )
        !
        !      CHARACTER     JOBVL, JOBVR
        !
        !      INTEGER       INFO, LDA, LDVL, LDVR, LWORK, N
        !
        !      DOUBLE        PRECISION RWORK( * )
        !
        !      COMPLEX*16    A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
        !                     W( * ), WORK( * )
        !
        ! PURPOSE
        !      ZGEEV computes for an N-by-N complex nonsymmetric matrix A,
        !      the eigenvalues and, optionally, the left and/or right
        !      eigenvectors.
        !
        !      The left eigenvectors of A are the same as the right eigen-
        !      vectors of A**H.  If u(j) and v(j) are the left and right
        !      eigenvectors, respectively, corresponding to the eigenvalue
        !      lambda(j), then (u(j)**H)*A = lambda(j)*(u(j)**H) and A*v(j)
        !      = lambda(j) * v(j).
        !
        !      The computed eigenvectors are normalized to have Euclidean
        !      norm equal to 1 and largest component real.
        !
        ! ARGUMENTS
        !      JOBVL   (input) CHARACTER*1
        !              = 'N': left eigenvectors of A are not computed;
        !              = 'V': left eigenvectors of are computed.
        !
        !      JOBVR   (input) CHARACTER*1
        !              = 'N': right eigenvectors of A are not computed;
        !              = 'V': right eigenvectors of A are computed.
        !
        !      N       (input) INTEGER
        !              The order of the matrix A. N >= 0.
        !
        !      A       (input/output) COMPLEX*16 array, dimension (LDA,N)
        !              On entry, the N-by-N matrix A.  On exit, A has been
        !              overwritten.
        !
        !      LDA     (input) INTEGER
        !              The leading dimension of the array A.  LDA >=
        !              max(1,N).
        !
        !      W       (output) COMPLEX*16 array, dimension (N)
        !              W contains the computed eigenvalues.
        !
        !      VL      (output) COMPLEX*16 array, dimension (LDVL,N)
        !              If JOBVL = 'V', the left eigenvectors u(j) are
        !              stored one after another in the columns of VL, in
        !              the same order as their eigenvalues.  If JOBVL =
        !              'N', VL is not referenced.  u(j) = VL(:,j), the j-th
        !              column of VL.
        !
        !      LDVL    (input) INTEGER
        !              The leading dimension of the array VL.  LDVL >= 1;
        !              if JOBVL = 'V', LDVL >= N.
        !
        !      VR      (output) COMPLEX*16 array, dimension (LDVR,N)
        !              If JOBVR = 'V', the right eigenvectors v(j) are
        !              stored one after another in the columns of VR, in
        !              the same order as their eigenvalues.  If JOBVR =
        !              'N', VR is not referenced.  v(j) = VR(:,j), the j-th
        !              column of VR.
        !
        !      LDVR    (input) INTEGER
        !              The leading dimension of the array VR.  LDVR >= 1;
        !              if JOBVR = 'V', LDVR >= N.
        !
        !      WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
        !              On exit, if INFO = 0, WORK(1) returns the optimal
        !              LWORK.
        !
        !      LWORK   (input) INTEGER
        !              The dimension of the array WORK.  LWORK >=
        !              max(1,2*N).  For good performance, LWORK must gen-
        !              erally be larger.
        !
        !      RWORK   (workspace) DOUBLE PRECISION array, dimension (2*N)
        !
        !      INFO    (output) INTEGER
        !              = 0:  successful exit
        !              < 0:  if INFO = -i, the i-th argument had an illegal
        !              value.
        !              > 0:  if INFO = i, the QR algorithm failed to com-
        !              pute all the eigenvalues, and no eigenvectors have
        !              been computed; elements and i+1:N of W contain
        !              eigenvalues which have converged.
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! use lapack95
        implicit none
        real*8,allocatable :: A(:,:)
        complex*16,allocatable :: RVEC(:,:),LAM(:,:),LVEC(:,:) !!! right and left eigvec
        integer*8 :: ind
        
        complex*16,allocatable :: Acpy(:,:)
        integer*8 :: N,LDA,LDVL,LDVR,LWORK
        integer*4 :: INFO
        character :: JOBVL = 'V', JOBVR = 'V'
        complex*16,allocatable :: W(:),VL(:,:),VR(:,:)
        complex*16,allocatable :: WORK(:),RWORK(:)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        N = size(A,dim=1);
        LDA = N; LDVL = N; LDVR = N; LWORK = 4*N;
        allocate(Acpy(N,N),W(N),VL(N,N),VR(N,N));
        
        Acpy = A + (0.0d0,0.0d0)
        
        allocate(WORK(4*N),RWORK(2*N));
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call ZGEEV(JOBVL, JOBVR, N, Acpy, LDA, W, VL, LDVL, VR, LDVR, &
                   WORK, LWORK, RWORK, INFO )
        
        if(INFO.NE.0)then
          STOP " Math_Eig ZGEEV Error "
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        LAM = 0.0d0
        DO ind = 1,N
           LAM(ind,ind) = (W(ind));
        END DO
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! ((VL)**(H))*A = LAM*((VL)**(H));
        !!! A*(VR) = LAM*(VR)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        LVEC = (VL);
        RVEC = (VR);
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        deallocate (Acpy,W,VL,VR);
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    END SUBROUTINE Math_Eig
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! LAPACK DGTRIF & DGETRI
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Math_Inv(A) RESULT(INV_A)
        !!! Compute the inverse matrix of real matrix A
        !!! Using DGETRF & DGETRI
        !!! See https://www.math.utah.edu/software/lapack/lapack-d/dgetri.html
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !NAME
        !  DGETRF - compute an LU factorization of a general M-by-N
        !  matrix A using partial pivoting with row interchanges
        !
        !SYNOPSIS
        !      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
        !
        !          INTEGER        INFO, LDA, M, N
        !
        !          INTEGER        IPIV( * )
        !
        !          DOUBLE         PRECISION A( LDA, * )
        !
        !PURPOSE
        !      DGETRF computes an LU factorization of a general M-by-N
        !      matrix A using partial pivoting with row interchanges.
        !
        !      The factorization has the form
        !         A = P * L * U
        !      where P is a permutation matrix, L is lower triangular with
        !      unit diagonal elements (lower trapezoidal if m > n), and U
        !      is upper triangular (upper trapezoidal if m < n).
        !
        !      This is the right-looking Level 3 BLAS version of the algo-
        !      rithm.
        !
        !ARGUMENTS
        !      M       (input) INTEGER
        !              The number of rows of the matrix A.  M >= 0.
        !
        !      N       (input) INTEGER
        !              The number of columns of the matrix A.  N >= 0.
        !
        !      A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
        !              On entry, the M-by-N matrix to be factored.  On
        !              exit, the factors L and U from the factorization A =
        !              P*L*U; the unit diagonal elements of L are not
        !              stored.
        !
        !      LDA     (input) INTEGER
        !              The leading dimension of the array A.  LDA >=
        !              max(1,M).
        !
        !      IPIV    (output) INTEGER array, dimension (min(M,N))
        !              The pivot indices; for 1 <= i <= min(M,N), row i of
        !              the matrix was interchanged with row IPIV(i).
        !
        !      INFO    (output) INTEGER
        !              = 0:  successful exit
        !              < 0:  if INFO = -i, the i-th argument had an illegal
        !              value
        !
        !              > 0:  if INFO = i, U(i,i) is exactly zero. The fac-
        !              torization has been completed, but the factor U is
        !              exactly singular, and division by zero will occur if
        !              it is used to solve a system of equations.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !NAME
        !DGETRI - compute the inverse of a matrix using the LU factorization computed by DGETRF
        !
        !SYNOPSIS
        !      SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
        !
        !      INTEGER        INFO, LDA, LWORK, N
        !
        !      INTEGER        IPIV( * )
        !
        !      DOUBLE         PRECISION A( LDA, * ), WORK( LWORK )
        !
        !PURPOSE
        !      DGETRI computes the inverse of a matrix using the LU factor-
        !      ization computed by DGETRF.
        !
        !      This method inverts U and then computes inv(A) by solving
        !      the system inv(A)*L = inv(U) for inv(A).
        !
        !ARGUMENTS
        !      N       (input) INTEGER
        !              The order of the matrix A.  N >= 0.
        !
        !      A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
        !              On entry, the factors L and U from the factorization
        !              A = P*L*U as computed by DGETRF.  On exit, if INFO =
        !              0, the inverse of the original matrix A.
        !
        !      LDA     (input) INTEGER
        !              The leading dimension of the array A.  LDA >=
        !              max(1,N).
        !
        !      IPIV    (input) INTEGER array, dimension (N)
        !              The pivot indices from DGETRF; for 1<=i<=N, row i of
        !              the matrix was interchanged with row IPIV(i).
        !
        !      WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
        !              On exit, if INFO=0, then WORK(1) returns the optimal
        !              LWORK.
        !
        !      LWORK   (input) INTEGER
        !              The dimension of the array WORK.  LWORK >= max(1,N).
        !              For optimal performance LWORK >= N*NB, where NB is
        !              the optimal blocksize returned by ILAENV.
        !
        !      INFO    (output) INTEGER
        !              = 0:  successful exit
        !              < 0:  if INFO = -i, the i-th argument had an illegal
        !              value
        !              > 0:  if INFO = i, U(i,i) is exactly zero; the
        !
        !              matrix is singular and its inverse could not be com-
        !              puted.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! use lapack95
        implicit none
        integer*4 :: INFO
        integer*8 :: LDA, LWORK, N 
        integer*4,allocatable :: IPIV(:)
        real*8,allocatable :: A(:,:), INV_A(:,:), Acpy(:,:),WORK(:)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        N = size(A,dim=1);
        LDA = N; LWORK = 4*N;
        
        allocate(IPIV(N));
        allocate(Acpy(N,N),INV_A(N,N),WORK(LWORK));
        
        Acpy = A;
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call DGETRF( N, N, Acpy, LDA, IPIV, INFO )
        
        if(INFO.NE.0)then
          STOP " Math_Inv DGETRF Error "
        end if
        
        call DGETRI( N, Acpy, LDA, IPIV, WORK, LWORK, INFO )
        
        if(INFO.NE.0)then
          STOP " Math_Inv DGETRI Error "
        end if
        
        INV_A = Acpy;
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        deallocate(IPIV);
        deallocate(Acpy,WORK);
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    END FUNCTION Math_Inv
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Math_Diagonal_ReEig_Mat(A,S,LAM,INV_S)
    !!!! diagnolize double precision matrix with real eigen value, be careful
    !!!! A = S * LAM * INV_S
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! use lapack95
        implicit none
        real*8,allocatable :: A(:,:),S(:,:),LAM(:,:),INV_S(:,:)
        complex*16,allocatable :: CS(:,:),CLAM(:,:)
        integer*8 :: N
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        N = size(A,dim=1);
        allocate(CS(N,N),CLAM(N,N)); CS = 0.0d0; CLAM = 0.0d0;
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call Math_Eig(A,CS,CLAM);
        S = dreal(CS);
        LAM = dreal(CLAM);
        INV_S = Math_Inv(S);
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        deallocate(CS,CLAM)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    END SUBROUTINE Math_Diagonal_ReEig_Mat
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Math_Diagonal_TwinReEig_Mat(A,S,LAM,INV_S)
    !!!! diagnolize double precision matrix with twin real eigen value, be careful
    !!!! e.g. +LAM(1), - LAM(1), +LAM(2), - LAM(2), ..., +LAM(N/2), -LAM(N/2);
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! A = S * LAM * INV_S
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! use lapack95
        implicit none
        real*8,allocatable :: A(:,:),S(:,:),LAM(:,:),INV_S(:,:)
        real*8,allocatable :: tmpS(:,:),tmpLAM(:,:),tmpLAMVec(:);
        complex*16,allocatable :: CS(:,:),CLAM(:,:)
        integer*8,allocatable :: TwinInd(:)
        integer*8 :: N,i
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        N = size(A,dim=1);
        allocate(tmpS(N,N),tmpLAM(N,N),tmpLAMVec(N),TwinInd(N));
        allocate(CS(N,N),CLAM(N,N)); 
        CS = 0.0d0; CLAM = 0.0d0;
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call Math_Eig(A,CS,CLAM);
        
        tmpS = dreal(CS);
        tmpLAM = dreal(CLAM);
        
        DO i =1,N
           tmpLAMVec(i) = tmpLAM(i,i); 
        END DO
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call Math_BubbleSort(tmpLAMVec,TwinInd,'T') !!! get the right twin ind
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! sort again
        LAM = 0.0d0;
        DO i = 1,N
           LAM(i,i) = tmpLAM(TwinInd(i),TwinInd(i));
        END DO
        
        S = 0.0d0;
        DO i = 1,N
           S(1:N,i) = tmpS(1:N,TwinInd(i));
        END DO
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        INV_S = Math_Inv(S);
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        deallocate(tmpS,tmpLAM,tmpLAMVec,TwinInd);
        deallocate(CS,CLAM); 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    END SUBROUTINE Math_Diagonal_TwinReEig_Mat
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Math_Interp2_3P
    !!! Two-Dim three-point interpolation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Math_InterP2_3P(Y,X,FD,YQ,XQ) RESULT(FDQ) 
    
    !!! LIKE MATLAB, Y as ROW, X as COL, be careful !!!!!!!!!
    
    !!! For ESLQ3, X = Y, Y = X, be careful
    !!! WARNING: Y(i) > Y(i+1) & X(j) < X(j+1)
    
        implicit none
        
        real*8,allocatable :: Y(:),X(:),FD(:,:),YQ(:),XQ(:),FDQ(:,:) !!! (Dim01,Dim02)
        
        real*8,allocatable :: OppoY(:) !!! opposite number of Y, to use ELSQ3      
        
        integer*8 :: YDim, YQ_Dim,XQ_Dim,I,J
        real*8 :: Val
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        YDim = size(Y)
        
        YQ_Dim = size(YQ)
        
        XQ_Dim = size(XQ)
        
        allocate(OppoY(YDim)); OppoY = (-1.0d0)*Y;
        
        allocate(FDQ(YQ_Dim,XQ_Dim))
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO I = 1, YQ_Dim
        
           DO J = 1, XQ_Dim
           
              CALL ESLQ3(OppoY,X,FD,(-1.0d0)*YQ(I),XQ(J),Val)
           
              FDQ(I,J) = Val
           
           END DO
           
        END DO
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CONTAINS
        
        SUBROUTINE ESLQ3(X,Y,Z,U,V,W)
            !!! mainly copyed from ShiLiang Xu (THU press,2nd Edition,1995)
            !!! it is a classcial numerical recipe in Chinese
            !!! see http://fcode.cn/resource_ebook-10-1.html
            !!! Section 5.16, Page 211
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            implicit none
        
            !!! Z = Z(Xi,Yi), WARNING : X(i) < X(i+1), Y(j) < Y(j+1)
        
	        !!! DIMENSION X(N),Y(M),Z(N,M),B(3)
            real*8,allocatable :: X(:),Y(:),Z(:,:)
            
	        !!! DOUBLE PRECISION X,Y,Z,U,V,W,B,HH
            real*8 :: B(3), U,V,W, HH
            
            integer*8 :: N,M,NN,MM,IP,IQ,I,J,K,L
            
            N = size(X); M = size(Y)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	        NN=3
	        IF (N.LE.3) THEN
	          IP=1
	          NN=N
	        ELSE IF (U.LE.X(2)) THEN
	          IP=1
	        ELSE IF (U.GE.X(N-1)) THEN
	          IP=N-2
	        ELSE
	          I=1
	          J=N
        10	  IF (IABS(I-J).NE.1) THEN
	            L=(I+J)/2
	            IF (U.LT.X(L)) THEN
	              J=L
	            ELSE
	              I=L
	            END IF
	            GOTO 10
	          END IF
	          IF (ABS(U-X(I)).LT.ABS(U-X(J))) THEN
	            IP=I-1
	          ELSE
	            IP=I
	          END IF
	        END IF
	        MM=3
	        IF (M.LE.3) THEN
	          IQ=1
	          MM=M
	        ELSE IF (V.LE.Y(2)) THEN
	          IQ=1
	        ELSE IF (V.GE.Y(M-1)) THEN
	          IQ=M-2
	        ELSE
	          I=1
	          J=M
        20	  IF (IABS(J-I).NE.1) THEN
	            L=(I+J)/2
	            IF (V.LT.Y(L)) THEN
	              J=L
	            ELSE
	              I=L
	            END IF
	            GOTO 20
	          END IF


	          IF (ABS(V-Y(I)).LT.ABS(V-Y(J))) THEN
	            IQ=I-1
	          ELSE
	            IQ=I
	          END IF
	        END IF
            
	        DO  I=1,NN
	          B(I)=0.0
	          DO  J=1,MM
	            HH=Z(IP+I-1,IQ+J-1)
	            DO  K=1,MM
	              IF (K.NE.J) THEN
	                HH=HH*(V-Y(IQ+K-1))/(Y(IQ+J-1)-Y(IQ+K-1))
	              END IF
        	    END DO
	            B(I)=B(I)+HH
        	  END DO
        	END DO
            
	        W=0.0
            
	        DO  I=1,NN
	          HH=B(I)
	          DO  J=1,NN
	            IF (J.NE.I) THEN
	              HH=HH*(U-X(IP+J-1))/(X(IP+I-1)-X(IP+J-1))
	            END IF
        	  END DO
	          W=W+HH
        	END DO
        
	        RETURN
            
	    END SUBROUTINE ESLQ3
    
    END FUNCTION Math_InterP2_3P
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Math_InterP1_2P
    !!! One-Dim two-point linear interpolation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Math_InterP1_2P(X,FD,XQ) RESULT(FDQ)
        !!! X(i) < X(i+1)
        implicit none
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer*8 :: Dim,I
        real*8,allocatable :: X(:),FD(:),XQ(:),FDQ(:)
        real*8 ::  Val
        
        Dim = size(XQ)
        allocate(FDQ(Dim))
        
        DO I = 1,Dim
        
           Val = Math_FindValue_In_MZSeq(X,FD,XQ(I),'A') 
           FDQ(I) = Val
        
        END DO
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    END FUNCTION Math_InterP1_2P
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Math_InterP1_3P
    !!! One-Dim three-point interpolation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Math_InterP1_3P(X,FD,XQ) RESULT (FDQ)
    
        implicit none
        integer*8 :: Dim,I
        real*8,allocatable :: X(:),XQ(:),FD(:),FDQ(:)
        real*8 :: Val
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Dim = size(XQ)
        allocate(FDQ(Dim));
        
        DO I = 1, Dim
        
           call ENLG3(X,FD,XQ(I),Val)
           
           FDQ(I) = Val
        
        END DO
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CONTAINS
        !!! mainly copyed from ShiLiang Xu (THU press,2nd Edition,1995)
        !!! it is a classcial numerical recipe in Chinese
        !!! see http://fcode.cn/resource_ebook-10-1.html
        !!! Section 5.3, Page 167
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE ENLG3(X,Y,T,Z)
            implicit none
            
            real*8,allocatable :: X(:),Y(:) !!! Y(X(N))
            integer*8 :: N,K,L,M,I,J
            real*8 :: T,Z,S
            
            N = size(X)
        
	        !!!DIMENSION X(N),Y(N)
	        !!!DOUBLE PRECISION X,Y,T,Z,S
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
	        Z=0.0
	        IF (N.LE.0) RETURN
	        IF (N.EQ.1) THEN
	          Z=Y(1)
	          RETURN
	        END IF
	        IF (N.EQ.2) THEN
	          Z=(Y(1)*(T-X(2))-Y(2)*(T-X(1)))/(X(1)-X(2))
	          RETURN
	        END IF
	        IF (T.LE.X(2)) THEN
	          K=1
	          M=3
	        ELSE IF (T.GE.X(N-1)) THEN
	          K=N-2
	          M=N
	        ELSE
	          K=1
	          M=N
        10	  IF (IABS(K-M).NE.1) THEN
	            L=(K+M)/2
	            IF (T.LT.X(L)) THEN
	              M=L
	            ELSE
	              K=L
	            END IF
	            GOTO 10
	          END IF
	          IF (ABS(T-X(K)).LT.ABS(T-X(M))) THEN
	            K=K-1
	          ELSE
	            M=M+1
	          END IF
	        END IF
            
	        Z=0.0
            
	        DO I=K,M
	          S=1.0
	          DO J=K,M
	            IF (J.NE.I) THEN
	              S=S*(T-X(J))/(X(I)-X(J))
	            END IF
        	  END DO
	          Z=Z+S*Y(I)
        	END DO
            
	        RETURN
            
        END SUBROUTINE ENLG3

    
    END FUNCTION Math_InterP1_3P
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Math_Smooth_5P
    !!!
    !!! cubical smoothing algorithm with five-point approximation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Math_Smooth_5P(X,M) RESULT(SX)
    
       implicit none
       
       real*8,allocatable :: X(:),SX(:)
       real*8,allocatable :: a(:)
       real*8,allocatable :: b(:)
       integer*8 :: m  !!! smooth times
       integer*8 :: n  !!! dim of A
       
       integer*8 :: j,k
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       n = size(X);
       allocate(SX(n),a(n),b(n))
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       a = X
       
       DO k = 1,m
          
          b(1) = ((69.0d0)*a(1) +(4.0d0)*(a(2) +a(4)) -(6.0d0)*a(3) -a(5)) /(70.0d0);
          
          b(2) = ((2.0d0)* (a(1) +a(5)) + (27.0d0)*a(2) + (12.0d0)*a(3) -(8.0d0)*a(4)) /(35.0d0);
        
          DO j =3,n-2
          
              b (j) = (-(3.0d0)*(a(j-2) +a(j+2)) +(12.0d0)*(a(j-1) +a(j+1)) +(17.0d0)*a(j))/(35.0d0);
              
          END DO
        
          b(n-1) = ((2.0d0)*(a(n) +a(n-4)) +(27.0d0)*a(n-1) +(12.0d0)*a(n-2) -(8.0d0)*a(n-3))/(35.0d0);
          
          b(n) = ((69.0d0)*a(n) + (4.0d0)* (a(n-1) +a(n-3)) -(6.0d0)*a(n-2) -a(n-4)) /(70.0d0);

          a = b
          
       END DO
       
       SX = b
    
    END FUNCTION Math_Smooth_5P
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Math_Smooth_9P
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Math_Smooth9P_2D(Mat,nsmooth)
    
        implicit none
        real*8,allocatable :: Mat(:,:) !!! row >=3, col>=3
        integer*8 :: nsmooth, IT, IY, IX, NY, NX
        real*8,allocatable :: HelpMat(:,:)
        real*8 :: Mwei=1.0d0, Bwei = 0.5d0, Cwei = 0.3d0
        real*8 :: TotalWei
        
        NY = size(Mat,1)-1
        NX = size(Mat,2)-1
        allocate(HelpMat(NY+1,NX+1)); HelpMat = 0.0d0;
        
        TotalWei = Mwei + (4.0d0)*Bwei + (4.0d0)*Cwei
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO IT= 1,nsmooth
           DO IY = 2,NY
              DO IX = 2,NX
              
                 HelpMat(IY,IX) = Mwei*Mat(IY,IX) + &
                     Bwei*( Mat(IY-1,IX)+ Mat(IY+1,IX) + Mat(IY,IX-1) + Mat(IY,IX+1) ) + &
                     Cwei*( Mat(IY-1,IX-1)+ Mat(IY+1,IX-1) + Mat(IY-1,IX+1) + Mat(IY+1,IX+1) )
                 
                 HelpMat(IY,IX) = HelpMat(IY,IX)/TotalWei
                 
    
              END DO
           END DO
           Mat = HelpMat;
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !!! set natural boundary condition 
           Mat(1,1:NX+1) = Mat(2,1:NX+1)
           Mat(NY+1,1:NX+1) = Mat(NY,1:NX+1)
        
           Mat(1:NY+1,1) = Mat(1:NY+1,2)
           Mat(1:NY+1,NX+1) = Mat(1:NY+1,NX)
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        END DO
        
        deallocate(HelpMat)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    END SUBROUTINE Math_Smooth9P_2D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Math_MeshGrid
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Math_MeshGrid(XX,YY,X,Y) 
    !!!! Just like meshgrid in Matlab
        real*8,allocatable :: X(:),XX(:,:)
        real*8,allocatable :: Y(:),YY(:,:)
        integer*8 :: NX,NY,IX,IY
        
        NX = size(X)
        NY = size(Y)
        
        DO IX = 1,NX
           YY(1:NY,IX) = Y
        END DO
        
        DO IY = 1,NY
           XX(IY,1:NX) = X
        END DO
    
    END SUBROUTINE Math_MeshGrid
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! UpWindXYAdv
    !!! LOW-ORDER UpWind Scheme
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION UpWindXYAdv(UL,UR,TL,TM,TR,DXY) RESULT(Adv)
    !!! UpWindXYAdv = [-u*(T)x = -(uT)x + T*(u)x] OR [v(T)y = -(vT)y + v*(T)y]
    !!! Staggered Arakawa-C Grid
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8 :: UL,UR
        real*8 :: TL,TM,TR
        real*8 :: DXY
        real*8 :: Adv
        
        !!! use low-order upwind scheme
        Adv = (-1.0d0/DXY) * (TM) * psgn(UR) + (-1.0d0/DXY) * (TR) * msgn(UR) + &
              (+1.0d0/DXY) * (TL) * psgn(UL) + (+1.0d0/DXY) * (TM) * msgn(UL)
              
        !!! be careful, this operation add non-conservative relationship into the model
        Adv = Adv + (1.0d0/DXY) * (TM) * (UR-UL)
    
    END FUNCTION UpWindXYAdv
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! UpWindXYFlux
    !!! LOW-ORDER UpWind Scheme
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION UpWindXYFlux(UL,UR,TL,TM,TR,DXY) RESULT(Flux)
    !!! UpWindXYFlux = [ -(uT)x ] OR [ -(vT)y ]
    !!! Staggered Arakawa-C Grid
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8 :: UL,UR
        real*8 :: TL,TM,TR
        real*8 :: Flux
        real*8 :: DXY
        
        !!! use low-order upwind scheme
        Flux = (-1.0d0/DXY) * (TM) * psgn(UR) + (-1.0d0/DXY) * (TR) * msgn(UR) + &
               (+1.0d0/DXY) * (TL) * psgn(UL) + (+1.0d0/DXY) * (TM) * msgn(UL)
        
    END FUNCTION UpWindXYFlux
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! SupBeeXYFlux
    !!!! SupBee TVD Scheme
    !!!! See Ocean Modelling for Beginners (Jochen Kampf,2009) P105 for more details
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION SupBeeXYFlux( UL, UR, TLL, TL, TM, TR, TRR, DXY, DT ) RESULT(Flux)
    
        real*8 :: UL,UR
        real*8 :: TLL,TL,TM,TR,TRR,DXY,DT
        real*8 :: Flux

        Flux = (-1.0d0/DXY) * psgn(UR) * SupBee(TL,TM,TR,TRR, CFLR = (DT/DXY)*abs(UR), SIG = +1.0d0)
        
        Flux = Flux + &
               (-1.0d0/DXY) * msgn(UR) * SupBee(TL,TM,TR,TRR, CFLR = (DT/DXY)*abs(UR), SIG = -1.0d0)
        
        Flux = Flux + &
               (+1.0d0/DXY) * psgn(UL) * SupBee(TLL,TL,TM,TR, CFLR = (DT/DXY)*abs(UL), SIG = +1.0d0)
        
        Flux = Flux + &
               (+1.0d0/DXY) * msgn(UL) * SupBee(TLL,TL,TM,TR, CFLR = (DT/DXY)*abs(UL), SIG = -1.0d0)
        
    
    END FUNCTION SupBeeXYFlux
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! SupBeeXYAdv
    !!!! SupBee TVD Scheme
    !!!! See Ocean Modelling for Beginners (Jochen Kampf,2009) P105 for more details
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION SupBeeXYAdv( UL, UR, TLL, TL, TM, TR, TRR, DXY, DT ) RESULT(Adv)
    
        real*8 :: UL,UR
        real*8 :: TLL,TL,TM,TR,TRR,DXY,DT
        real*8 :: Adv

        Adv = (-1.0d0/DXY) * psgn(UR) * SupBee(TL,TM,TR,TRR, CFLR = (DT/DXY)*abs(UR), SIG = +1.0d0)
        
        Adv = Adv + &
               (-1.0d0/DXY) * msgn(UR) * SupBee(TL,TM,TR,TRR, CFLR = (DT/DXY)*abs(UR), SIG = -1.0d0)
        
        Adv = Adv + &
               (+1.0d0/DXY) * psgn(UL) * SupBee(TLL,TL,TM,TR, CFLR = (DT/DXY)*abs(UL), SIG = +1.0d0)
    
        Adv = Adv + &
               (+1.0d0/DXY) * msgn(UL) * SupBee(TLL,TL,TM,TR, CFLR = (DT/DXY)*abs(UL), SIG = -1.0d0)
        
        !!! be careful, this operation add non-conservative relationship into the model
        Adv = Adv + (1.0d0/DXY) * (TM) * (UR-UL)
    
    END FUNCTION SupBeeXYAdv
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! SupBee
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION SupBee( BLL,BL,BR,BRR,CFLR,SIG ) RESULT(SupB)
    
        real*8 :: BLL,BL,BR,BRR
        real*8 :: CFLR !!! CFL number = (DT/DX)*abs(VLR), VLR: velocity at the rightside interface
        real*8 :: SupB
        real*8 :: SIG !!! +(1.0d0) or -(1.0d0)
        
        real*8 :: Rplus,Rminus
        real*8 :: epi = 1.0d-14; !!! to avoid overflow of denominator
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        Rplus  = (BL-BLL)/(BR-BL + epi)
        Rminus = (BRR-BR)/(BR-BL + epi)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        SupB = &
        (+psgn(SIG)) * (BL + 0.5d0 * Bee_Limit(Rplus) * (1.0d0-CFLR) * (BR-BL) ) + &
        (-msgn(SIG)) * (BR - 0.5d0 * Bee_Limit(Rminus)* (1.0d0-CFLR) * (BR-BL) );
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    END FUNCTION SupBee
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! TVD limiter 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Bee_Limit(R)RESULT(PHI)
    
       real*8 :: R,PHI
       
       !!! unsmooth 
       PHI = max(0.0d0,min(2.0d0*R,1.0d0),min(R,2.0d0))
       
       !!! smooth verision (slow for compuation)
       !!! PHI = psgn( smmax( UbdCtrl(2.0d0*R,1.0d0), UbdCtrl(R,2.0d0) ) );
        
    END FUNCTION Bee_Limit
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
END MODULE MO_MATHTOOLS
