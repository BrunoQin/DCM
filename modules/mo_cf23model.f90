MODULE MO_CF23Model
    
USE MO_NCTOOLS
USE MO_MATHTOOLS
USE MO_TYPEDEF

IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 
!!! The Fortran Version of 
!!! Nan Chen & Fang Xianghui (2023) (JAMES) 
!!! A Simple Multiscale Intermediate Coupled Stochastic Model for El Nino Diversity and Complexity
!!! https://doi.org/10.1029/2022MS003469
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CFModel_Initialization (CFM_Par, CFM_Stat,TimTick)

    TYPE(CFModel_Param) :: CFM_Par
    TYPE(CFModel_Stat)  :: CFM_Stat
    TYPE(TimeTick) :: TimTick  
    
    integer*8 :: NAtm,NOcn,AI_NX,AI_NY
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (TimTick.Run_Cycle_ID == 1 ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    NAtm = CFM_Par.NAtm
    NOcn = CFM_Par.NOcn
    AI_NX = CFM_Par.AI_NX
    AI_NY = CFM_Par.AI_NY
    
    allocate(CFM_Par.OLon(NOcn));allocate(CFM_Par.OLat(AI_NY));
    allocate(CFM_Par.OReX(NOcn));allocate(CFM_Par.OReY(AI_NY));
    
    allocate(CFM_Par.ALon(NAtm));allocate(CFM_Par.ALat(AI_NY));
    allocate(CFM_Par.AReX(NAtm));allocate(CFM_Par.AReY(AI_NY));
    
    !!! ocn,atm lon
    !!! CFM_Par.OLon = Math_Linspace(117.0d0,282.0d0,NOcn);
    CFM_Par.OLon = Math_Linspace(CFM_Par.OWBD,CFM_Par.OEBD,NOcn);
    CFM_Par.ALon = CFM_Par.OWBD + Math_Linspace(0.0d0,360.0d0-(360.0d0/dble(NAtm)),NAtm);
    
    !!! for (288x64) AI Model + ELBOM
    CFM_Par.OLat = Math_Linspace(CFM_Par.OSBD,CFM_Par.ONBD,AI_NY);
    CFM_Par.ALat = Math_Linspace(CFM_Par.OSBD,CFM_Par.ONBD,AI_NY);
    
    CFM_Par.OReX = (OneDegLen) * CFM_Par.OLon;
    CFM_Par.OReY = (OneDegLen) * CFM_Par.OLat;
    CFM_Par.AReX = (OneDegLen) * CFM_Par.ALon;
    CFM_Par.AReY = (OneDegLen) * CFM_Par.ALat;
    
    allocate(CFM_Par.ONdReY(AI_NY)); 
    allocate(CFM_Par.ANdReY(AI_NY));
    
    CFM_Par.ONdReY = (CFM_Par.OReY) / sqrt(CFM_Par.COcn / EarthBeta);
    CFM_Par.ANdReY = (CFM_Par.AReY) / sqrt(CFM_Par.CAtm / EarthBeta);
    
    !!! get meridional modes
    !!! phi_0 = 1/pi^(1/4) * exp(-yy.^2/2);
    !!! phi_2 = 1/pi^(1/4)/sqrt(8) * (4 * yy.^2 - 2) .* exp(-yy.^2/2);
    allocate (CFM_Par.OPHI_0(AI_NY));
    allocate (CFM_Par.OPHI_2(AI_NY));
    allocate (CFM_Par.APHI_0(AI_NY));
    allocate (CFM_Par.APHI_2(AI_NY));
    
    CFM_Par.OPHI_0 = (1.0d0 /sqrt(sqrt(pi)) ) * exp( -(0.5d0) * (CFM_Par.ONdReY)**(2.0d0) );
    CFM_Par.APHI_0 = (1.0d0 /sqrt(sqrt(pi)) ) * exp( -(0.5d0) * (CFM_Par.ANdReY)**(2.0d0) );
    
    CFM_Par.OPHI_2 = (1.0d0 /( sqrt(sqrt(pi)) * sqrt(8.0d0) ) ) * & 
                     (4.0d0 * (CFM_Par.ONdReY)**(2.0d0) - 2.0d0 ) * &
                     exp( -(0.5d0) * (CFM_Par.ONdReY)**(2.0d0) );
    
    CFM_Par.APHI_2 = (1.0d0 /( sqrt(sqrt(pi)) * sqrt(8.0d0) ) ) * & 
                     (4.0d0 * (CFM_Par.ANdReY)**(2.0d0) - 2.0d0 ) * &
                     exp( -(0.5d0) * (CFM_Par.ANdReY)**(2.0d0) );
    
    
    CFM_Par.Eq_OPHI_0 = (1.0d0 /sqrt(sqrt(pi)) );
    CFM_Par.Eq_APHI_0 = (1.0d0 /sqrt(sqrt(pi)) );
    
    CFM_Par.Eq_OPHI_2 = (1.0d0 /( sqrt(sqrt(pi)) * sqrt(8.0d0) ) ) * ( - 2.0d0 );
    CFM_Par.Eq_APHI_2 = (1.0d0 /( sqrt(sqrt(pi)) * sqrt(8.0d0) ) ) * ( - 2.0d0 );
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    call NC_Read_Mul_IT(Grid_Path,CFModel_Param_FName,'Ka_temp',CFM_Par.Ka_temp,0,0,'A');
    call NC_Read_Mul_IT(Grid_Path,CFModel_Param_FName,'Ra_temp',CFM_Par.Ra_temp,0,0,'A');
    
    call NC_Read_Mul_IT(Grid_Path,CFModel_Param_FName,'RandN',CFM_Par.RandN,0,0,'A');
    call NC_Read_Mul_IT(Grid_Path,CFModel_Param_FName,'Rd_Save',CFM_Par.Rd_Save,0,0,'A');
    
    call NC_Read_One_IT(Grid_Path,CFModel_Param_FName,'m',CFM_Par.IdM,1,'A');
    call NC_Read_One_IT(Grid_Path,CFModel_Param_FName,'lambda',CFM_Par.IdLambda,1,'A');
    call NC_Read_Mul_IT(Grid_Path,CFModel_Param_FName,'sgm',CFM_Par.IdSgm,0,0,'A');
    
    call NC_Read_Mul_IT(Grid_Path,CFModel_Param_FName,'Lin_M',CFM_Par.Lin_M,0,0,'A');
    call NC_Read_Mul_IT(Grid_Path,CFModel_Param_FName,'eta',CFM_Par.eta,0,0,'A');
    call NC_Read_Mul_IT(Grid_Path,CFModel_Param_FName,'eta2',CFM_Par.eta2,0,0,'A');
    call NC_Read_Mul_IT(Grid_Path,CFModel_Param_FName,'sp',CFM_Par.sp,0,0,'A');
        

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! compute global_sp
    allocate(CFM_Par.glob_sp(NAtm)); CFM_Par.glob_sp = 0.0d0;
    
    CFM_Par.glob_sp = 0.0d0;
    
    CFM_Par.glob_sp = &
    exp(-(45.0d0)*( ( CFM_Par.dx * Math_Linspace(1.0d0,128.0d0,128) - (CFM_Par.Lo)/4.0d0) )**(2.0d0) );
    
    CFM_Par.glob_sp (NAtm-30:NAtm) =  &
    exp(-(45.0d0)*( ( CFM_Par.dx * Math_Linspace(-30.0d0,0.0d0,31) - (CFM_Par.Lo)/4.0d0) )**(2.0d0) );
    
    
    !!! write(*,*),'CFM_Par.IdM =',CFM_Par.IdM
    !!! change CFM_Par.IdM
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! Initialize CFM_Stat
    allocate(CFM_Stat.KAtm(NAtm)); CFM_Stat.KAtm = 0.0d0;
    allocate(CFM_Stat.RAtm(NAtm)); CFM_Stat.RAtm = 0.0d0;
    
    allocate(CFM_Stat.KOcn(NOcn));     CFM_Stat.KOcn = 0.0d0;
    allocate(CFM_Stat.ROcn(NOcn));     CFM_Stat.ROcn = 0.0d0;
    allocate(CFM_Stat.Nd_Ocn_T(NOcn)); CFM_Stat.Nd_Ocn_T = 0.0d0;
    
    !!! assemble state Vector
    !!! Note
    !!! CFM_Stat.UVec(1:NOcn) = CFM_Stat.KOcn
    !!! CFM_Stat.UVec(NOcn+1 : 2*NOcn) = CFM_Stat.ROcn
    !!! CFM_Stat.Ocn_T = OPHI_0_Eq * Dim_Ocn_T * UVec(2*NOcn+1 : 3*NOcn)
    allocate(CFM_Stat.UVec(3*NOcn,1)); CFM_Stat.UVec = 0.0d0;
    
    
    !!! dimensional variable
    allocate(CFM_Stat.Ocn_U(NOcn)); CFM_Stat.Ocn_U = 0.0d0;
    allocate(CFM_Stat.Ocn_H(NOcn)); CFM_Stat.Ocn_H = 0.0d0;
    allocate(CFM_Stat.Ocn_T(NOcn)); CFM_Stat.Ocn_T = 0.0d0;
    allocate(CFM_Stat.Atm_U(NOcn)); CFM_Stat.Atm_U = 0.0d0;
    
    allocate(CFM_Stat.Ocn_UU(NOcn,AI_NY)); CFM_Stat.Ocn_UU = 0.0d0;
    allocate(CFM_Stat.Ocn_HH(NOcn,AI_NY)); CFM_Stat.Ocn_HH = 0.0d0;
    allocate(CFM_Stat.Ocn_TT(NOcn,AI_NY)); CFM_Stat.Ocn_TT = 0.0d0;
    allocate(CFM_Stat.Atm_UU(NOcn,AI_NY)); CFM_Stat.Atm_UU = 0.0d0;
    
    allocate(CFM_Stat.EvoMat(NOcn*3,NOcn*3)); CFM_Stat.EvoMat = 0.0d0;
    allocate(CFM_Stat.FrcVec(NOcn*3,1));      CFM_Stat.FrcVec = 0.0d0;
    allocate(CFM_Stat.Ex_FrcVec(NOcn*3,1));   CFM_Stat.Ex_FrcVec = 0.0d0;
    
    allocate(CFM_Stat.Glob_Atm_U(NAtm)); CFM_Stat.Glob_Atm_U = 0.0d0;
    allocate(CFM_Stat.Glob_WWB_U(NAtm)); CFM_Stat.Glob_WWB_U = 0.0d0;
    
    allocate(CFM_Stat.Ex_Ocn_Str(NOcn)); CFM_Stat.Ex_Ocn_Str = 0.0d0;
    
    !!!! simplified heat budget of CF23 (unit: dC/s not dC/dt*, be careful)
    allocate(CFM_Stat.TEND(NOcn)); CFM_Stat.TEND = 0.0d0;
    allocate(CFM_Stat.ZAFK(NOcn)); CFM_Stat.ZAFK = 0.0d0;
    allocate(CFM_Stat.THFK(NOcn)); CFM_Stat.THFK = 0.0d0;
    allocate(CFM_Stat.RESQ(NOcn)); CFM_Stat.RESQ = 0.0d0;
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!! monthly mean state of anomalies
    allocate(CFM_Stat.MonMe_KOcn(NOcn)); CFM_Stat.MonMe_KOcn = 0.0d0;
    allocate(CFM_Stat.MonMe_ROcn(NOcn)); CFM_Stat.MonMe_ROcn = 0.0d0;
    allocate(CFM_Stat.MonMe_Nd_Ocn_T(Nocn)); CFM_Stat.MonMe_Nd_Ocn_T = 0.0d0;
    
    allocate(CFM_Stat.MonMe_Ocn_U(NOcn)); CFM_Stat.MonMe_Ocn_U = 0.0d0;
    allocate(CFM_Stat.MonMe_Ocn_H(NOcn)); CFM_Stat.MonMe_Ocn_H = 0.0d0;
    allocate(CFM_Stat.MonMe_Ocn_T(NOcn)); CFM_Stat.MonMe_Ocn_T = 0.0d0;
    allocate(CFM_Stat.MonMe_Atm_U(NOcn)); CFM_Stat.MonMe_Atm_U = 0.0d0;
    
    allocate(CFM_Stat.MonMe_TEND(NOcn)); CFM_Stat.MonMe_TEND = 0.0d0;
    allocate(CFM_Stat.MonMe_ZAFK(NOcn)); CFM_Stat.MonMe_ZAFK = 0.0d0;
    allocate(CFM_Stat.MonMe_THFK(NOcn)); CFM_Stat.MonMe_THFK = 0.0d0;
    allocate(CFM_Stat.MonMe_RESQ(NOcn)); CFM_Stat.MonMe_RESQ = 0.0d0;
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Assigan initial value of state
    
    !!! default initial value of interdecadal variability & intraseasonal wind burst
    CFM_Stat.IntDec = 0.5d0 
    CFM_Stat.ItrAmp = 0.0d0

    
    !!! KOcn, ROcn, Nd_Ocn_T
    !!! MATLAB code:
    !!! tmp = sin([1/No:1/No:1]*2*pi);
    !!! u0(0*No+1:3*No,:) = 2*[tmp,-tmp/2,tmp];
    CFM_Stat.KOcn     = (2.0d0)  * sin ( Math_Linspace( 1.0d0/dble(NOcn), 1.0d0, NOcn ) * (2.0d0*pi) );
    CFM_Stat.ROcn     = (-0.5d0) * CFM_Stat.KOcn;
    CFM_Stat.Nd_Ocn_T = CFM_Stat.KOcn;
    
    CFM_Stat.UVec(         1:   NOcn       , 1) = CFM_Stat.KOcn;
    CFM_Stat.UVec(  NOcn + 1:   NOcn + NOcn, 1) = CFM_Stat.ROcn;
    CFM_Stat.UVec(2*NOcn + 1: 2*NOcn + NOcn, 1) = CFM_Stat.Nd_Ocn_T;
    
    CFM_Stat.EvoMat = CFM_Par.Lin_M;
    CFM_Stat.FrcVec = 0.0d0;
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! arbitrary intialization
    !!! CFM_Stat.UVec(         1: 2*NOcn       , 1) = (0.0d0)
    !!! CFM_Stat.UVec(2*NOcn + 1: 2*NOcn + NOcn, 1) = (-0.5d0)/(CFM_Par.Dim_Ocn_T);
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! T4 = mean(u0(T4_node));
    CFM_Stat.T4 = sum(  CFM_Stat.UVec(2*NOcn + CFM_Par.T4_WInd : 2*NOcn + CFM_Par.T4_EInd, 1))/&
                  dble( CFM_Par.T4_EInd - CFM_Par.T4_WInd + 1);
    
    !!! write(*,*),'CFM_Stat.T4 = ',CFM_Stat.T4
    
    !!! write(*,*),'CFM_Par.Re2Nd_Str = ',CFM_Par.Re2Nd_Str
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*),'------------------------------------------------------------------------------'
    write(*,*),'Chen-Fang 2023 Intermeidate Coupled Model Initialized'
    write(*,*),'------------------------------------------------------------------------------'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END IF !!! TimTick.Run_Cycle_ID == 1
    

END SUBROUTINE CFModel_Initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CFModel_Write_ReStart(ResPath,CFM_Par,CFM_Stat,TimTick);

     character(len=*) :: ResPath
     TYPE(CFModel_Param) :: CFM_Par
     TYPE(CFModel_Stat)  :: CFM_Stat
     TYPE(TimeTick) :: TimTick  
    
     
     character(len=300) :: CFM_ResName = 'CF23_ReStart.nc';
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     real*8,allocatable :: res_time_tick(:) !!! tmp timetick size == 1
     character(len=100) :: time_unit   !!! time_unit
     
     allocate(res_time_tick(1)); res_time_tick(1) = dble(TimTick.Run_AbsIDay);
     
     !!! to prevent misjudge of CDO (sub half-day)
     res_time_tick = res_time_tick - 0.5d0
     
     time_unit = &
     'days since '//trim(Math_YMD_To_Date(Start_CYear,Start_CMonth,Start_CDay))//' 00:00:0.0';
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     write(*,*),'------------------------------------------------------------------------------'
     write(*,*),'Writing Chen-Fang 2023 Model ReStart Files ...'
     
     call  NC_Create_File( ResPath,CFM_ResName, &
                           CFM_Par.OLon, CFM_Par.OLat, res_time_tick, &
                           'lon', 'lat','time', time_unit, 'X', Leap_Calendar_OnOff)
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     !!! be careful, in restart file, numerical precision should be double
     call NC_Def_Var(ResPath,CFM_ResName,'KOcn','none','lon','time',Double_FLAG=1);
     call NC_Def_Var(ResPath,CFM_ResName,'ROcn','none','lon','time',Double_FLAG=1);
     call NC_Def_Var(ResPath,CFM_ResName,'Nd_Ocn_T','none','lon','time',Double_FLAG=1);
     call NC_Def_Var(ResPath,CFM_ResName,'IntDec','none','time',Double_FLAG=1);
     call NC_Def_Var(ResPath,CFM_ResName,'ItrAmp','none','time',Double_FLAG=1);
     
     call NC_Def_Var(ResPath,CFM_ResName,'Res_AbsIT','time step','time',Double_FLAG=1);
     call NC_Def_Var(ResPath,CFM_ResName,'CYear','year','time',Double_FLAG=1);
     call NC_Def_Var(ResPath,CFM_ResName,'CMonth','month','time',Double_FLAG=1);
     call NC_Def_Var(ResPath,CFM_ResName,'CDay','day','time',Double_FLAG=1);
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     call NC_Write_One_IT(ResPath,CFM_ResName,'KOcn',CFM_Stat.KOcn,1);
     call NC_Write_One_IT(ResPath,CFM_ResName,'ROcn',CFM_Stat.ROcn,1);
     call NC_Write_One_IT(ResPath,CFM_ResName,'Nd_Ocn_T',CFM_Stat.Nd_Ocn_T,1);
     call NC_Write_One_IT(ResPath,CFM_ResName,'IntDec',CFM_Stat.IntDec,1);
     call NC_Write_One_IT(ResPath,CFM_ResName,'ItrAmp',CFM_Stat.ItrAmp,1);
     
     !!! BE CAREFUL : TimeTick.Run_AbsIT >> Res_AbsIT
     call NC_Write_One_IT(ResPath,CFM_ResName,'Res_AbsIT',dble(TimTick.Run_AbsIT),1);
     call NC_Write_One_IT(ResPath,CFM_ResName,'CYear',dble(TimTick.CYear),1);
     call NC_Write_One_IT(ResPath,CFM_ResName,'CMonth',dble(TimTick.CMonth),1);
     call NC_Write_One_IT(ResPath,CFM_ResName,'CDay',dble(TimTick.CDay),1);
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     !!! make sure not empty file output 
     write(*,*),'Eq_T4       = '//trim(Num2Str(CFM_Stat.Eq_T4))//&
                '  Eq_T34       = '//trim(Num2Str(CFM_Stat.Eq_T34))//&
                '  Eq_T3        = '//trim(Num2Str(CFM_Stat.Eq_T3));
     write(*,*),''
     
     write(*,*),'Chen-Fang 2023 Model ReStart Files Writing Completed'
     write(*,*),'------------------------------------------------------------------------------'
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
END SUBROUTINE CFModel_Write_ReStart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CFModel_Read_ReStart(ResPath,CFM_Par,CFM_Stat,TimTick);

     character(len=*) :: ResPath
     TYPE(CFModel_Param) :: CFM_Par
     TYPE(CFModel_Stat)  :: CFM_Stat
     TYPE(TimeTick) :: TimTick  
     
     
     real*8 :: db_Res_AbsIT, db_CYear, db_CMonth, db_CDay
     !!! double precision absolute integration IT from one Run for restart
     character(len=300) :: CFM_ResName = 'CF23_ReStart.nc';
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     write(*,*)'------------------------------------------------------------------------------'
     write(*,*),'Chen-Fang 2023 Model Reading ReStart...'
     
     call NC_Read_One_IT(ResPath,CFM_ResName,'KOcn',CFM_Stat.KOcn,1,'N');
     call NC_Read_One_IT(ResPath,CFM_ResName,'ROcn',CFM_Stat.ROcn,1,'N');
     call NC_Read_One_IT(ResPath,CFM_ResName,'Nd_Ocn_T',CFM_Stat.Nd_Ocn_T,1,'N');
     call NC_Read_One_IT(ResPath,CFM_ResName,'IntDec',CFM_Stat.IntDec,1,'N');
     call NC_Read_One_IT(ResPath,CFM_ResName,'ItrAmp',CFM_Stat.ItrAmp,1,'N');
     
     !!! assemble UVec
     CFM_Stat.UVec(               1 :   CFM_Par.NOcn,1) = CFM_Stat.KOcn;
     CFM_Stat.UVec(  CFM_Par.NOcn+1 : 2*CFM_Par.NOcn,1) = CFM_Stat.ROcn;
     CFM_Stat.UVec(2*CFM_Par.NOcn+1 : 3*CFM_Par.NOcn,1) = CFM_Stat.Nd_Ocn_T;
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!! The following operations can be completed by ELBOM
     IF( CFModel_RunAlone_OnOff == 1 ) THEN
     
         call NC_Read_One_IT(ResPath,CFM_ResName,'Res_AbsIT',db_Res_AbsIT,1,'N');
         call NC_Read_One_IT(ResPath,CFM_ResName,'CYear',db_CYear,1,'N');
         call NC_Read_One_IT(ResPath,CFM_ResName,'CMonth',db_CMonth,1,'N');
         call NC_Read_One_IT(ResPath,CFM_ResName,'CDay',db_CDay,1,'N');
     
         TimTick.Res_AbsIT = int8(db_Res_AbsIT);
         TimTick.CYear     = int8(db_CYear);
         TimTick.CMonth    = int8(db_CMonth);
         TimTick.CDay      = int8(db_CDay);
     
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!! Clip_Run, directly assign days of run 
         IF ( Clip_Run_OnOff == 1 ) THEN
     
            TimTick.Run_AbsNT   = Clip_Run_NDay * (Day_NT) + TimTick.Res_AbsIT
            TimTick.Run_AbsNDay = Clip_Run_NDay + ceiling(dble(TimTick.Res_AbsIT)/dble(Day_NT))
     
         END IF
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
         write(6,120) TimTick.Res_AbsIT
         120  format(' TimTick.Res_AbsIT   = ',I10)
         write(6,130) TimTick.Run_AbsNT
         130  format(' TimTick.Res_AbsNT   = ',I10)
         write(6,140) TimTick.Run_AbsNDay
         140  format(' TimTick.Res_AbsNDay = ',I10)
     
         write(6,121) TimTick.CYear
         121  format(' TimTick.CYear       = ',I8)
         write(6,122) TimTick.CMonth
         122  format(' TimTick.CMonth      = ',I8)
         write(6,123) TimTick.CDay
         123  format(' TimTick.CDay        = ',I8)
     
     ENDIF !!! CFModel_RunAlone_OnOff
     
     write(*,*),'Chen-Fang 2023 Model ReStart Reading Completed'
     write(*,*),'---------------------------------------------------------------------------'
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
END SUBROUTINE CFModel_Read_ReStart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CFModel_ReInit_Stat(CFM_Par,CFM_Stat,TimTick)

     TYPE(CFModel_Param) :: CFM_Par
     TYPE(CFModel_Stat)  :: CFM_Stat
     TYPE(TimeTick) :: TimTick  
     
     integer*8 :: NOcn
     
     NOcn = CFM_Par.NOcn
     
     
     !!! default initial value of interdecadal variability & intraseasonal wind burst
     CFM_Stat.IntDec = 0.5d0 
     CFM_Stat.ItrAmp = 0.0d0

    
     !!! KOcn, ROcn, Nd_Ocn_T
     !!! MATLAB code:
     !!! tmp = sin([1/No:1/No:1]*2*pi);
     !!! u0(0*No+1:3*No,:) = 2*[tmp,-tmp/2,tmp];
     CFM_Stat.KOcn     = (2.0d0)  * sin ( Math_Linspace( 1.0d0/dble(NOcn), 1.0d0, NOcn ) * (2.0d0*pi) );
     CFM_Stat.ROcn     = (-0.5d0) * CFM_Stat.KOcn;
     CFM_Stat.Nd_Ocn_T = CFM_Stat.KOcn;
    
     CFM_Stat.UVec(         1:   NOcn       , 1) = CFM_Stat.KOcn;
     CFM_Stat.UVec(  NOcn + 1:   NOcn + NOcn, 1) = CFM_Stat.ROcn;
     CFM_Stat.UVec(2*NOcn + 1: 2*NOcn + NOcn, 1) = CFM_Stat.Nd_Ocn_T;
     

END SUBROUTINE CFModel_ReInit_Stat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CFModel_Create_OutPut(OutPath,CFM_Par,CFM_Stat,TimTick);

    character(len=*) :: OutPath
    TYPE(CFModel_Param) :: CFM_Par
    TYPE(CFModel_Stat)  :: CFM_Stat
    TYPE(TimeTick) :: TimTick  
    
    character(len=300) :: CFM_OutName
    character(len=300) :: CFM_Atm_OutName
    
    integer*8 :: IT
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
    real*8,allocatable :: out_time_tick(:) !!! tmp timetick size == OutPut_N_Event
    character(len=100) :: time_unit   !!! time_unit
     
    allocate(out_time_tick(OutPut_N_Event))
     
    !!! value for output_time_tick
    IF(TimTick.Run_AbsIDay == 0 ) THEN
    DO  IT = 1, OutPut_N_Event
        out_time_tick (IT) =  PutData_N_Day * IT 
    END DO
    ELSE
    !!! because we have entering new day for TimTick
    DO  IT = 1, OutPut_N_Event
        out_time_tick (IT) =  PutData_N_Day * IT + (TimTick.Run_AbsIDay - 1)
    END DO
    END IF
    
    
    !!! to prevent misjudge of CDO (sub half-day)
    out_time_tick = out_time_tick - 0.5d0
     
    !!! time_unit
    time_unit = &
    'days since '//trim(Math_YMD_To_Date(Start_CYear,Start_CMonth,Start_CDay))//' 00:00:0.0';
     
    !!! CFM_OutName
    CFM_OutName = 'CF23_'//trim(TimTick.Out_Date)//'_'//trim(Num2Str(OutPut_Day_Span))//'day_'// &
                  trim(Num2Str(PutData_N_Day))//'day_R'//trim(Num2Str(TimTick.Run_Cycle_ID))//'.nc'
    
    IF(Clip_Run_OnOff == 0)THEN
    
    CFM_Atm_OutName = 'CF23_Atm_'//trim(TimTick.Out_Date)//'_'//trim(Num2Str(OutPut_Day_Span))//'day_'// &
                       trim(Num2Str(PutData_N_Day))//'day_R'//trim(Num2Str(TimTick.Run_Cycle_ID))//'.nc'
    
    END IF !!! IF(Clip_Run_OnOff == 0)
     
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*),'----------------------------------------------------------------------------'
    write(*,*),'Creating Chen-Fang 2023 Model OutPut Files ...'

    !!!! Ocn Field of CF23
    
    call NC_Create_File( OutPath,CFM_OutName, &
                         CFM_Par.OLon, CFM_Par.OLat,out_time_tick, &
                         'lon', 'lat','time', time_unit, 'X', Leap_Calendar_OnOff)
     
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
    !!! OutPut timetick information
    call NC_Def_Var(OutPath,CFM_OutName,'Res_AbsIT','time step','time',Double_FLAG=1);
    call NC_Def_Var(OutPath,CFM_OutName,'CYear','year','time',Double_FLAG=1);
    call NC_Def_Var(OutPath,CFM_OutName,'CMonth','month','time',Double_FLAG=1);
    call NC_Def_Var(OutPath,CFM_OutName,'CDay','day','time',Double_FLAG=1);
     
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! non-dimensional variable
    call NC_Def_Var(OutPath,CFM_OutName,'KOcn','none','lon','time',Double_FLAG=1)
    call NC_Def_Var(OutPath,CFM_OutName,'ROcn','none','lon','time',Double_FLAG=1)
    call NC_Def_Var(OutPath,CFM_OutName,'Nd_Ocn_T','none','lon','time',Double_FLAG=1)
    call NC_Def_Var(OutPath,CFM_OutName,'IntDec','none','time',Double_FLAG=1);
    call NC_Def_Var(OutPath,CFM_OutName,'ItrAmp','none','time',Double_FLAG=1);
    
    !!! dimensional variable 
    call NC_Def_Var(OutPath,CFM_OutName,'Ocn_U','m/s','lon','time',Double_FLAG=1)
    call NC_Def_Var(OutPath,CFM_OutName,'Ocn_H','m','lon','time',Double_FLAG=1)
    call NC_Def_Var(OutPath,CFM_OutName,'Ocn_T','dC','lon','time',Double_FLAG=1)
    call NC_Def_Var(OutPath,CFM_OutName,'Atm_U','m/s','lon','time',Double_FLAG=1)
    call NC_Def_Var(OutPath,CFM_OutName,'Ex_Ocn_Str','N/m^2','lon','time',Double_FLAG=1);
    
    call NC_Def_Var(OutPath,CFM_OutName,'ZAFK','dC/s','lon','time',Double_FLAG=1)
    call NC_Def_Var(OutPath,CFM_OutName,'THFK','dC/s','lon','time',Double_FLAG=1)
    call NC_Def_Var(OutPath,CFM_OutName,'RESQ','dC/s','lon','time',Double_FLAG=1)
     
    IF( CFModel_OutPut_2D_OnOff == 1) THEN
    
        call NC_Def_Var(OutPath,CFM_OutName,'Ocn_UU','m/s','lon','lat','time',Double_FLAG=1)
        call NC_Def_Var(OutPath,CFM_OutName,'Ocn_HH','m','lon','lat','time',Double_FLAG=1)
        call NC_Def_Var(OutPath,CFM_OutName,'Ocn_TT','dC','lon','lat','time',Double_FLAG=1)
        call NC_Def_Var(OutPath,CFM_OutName,'Atm_UU','m/s','lon','lat','time',Double_FLAG=1)
    
    END IF !!! CFModel_OutPut_2D_OnOff
    
    !!! Nino T4,T34,T3 value on the equator
    call NC_Def_Var(OutPath,CFM_OutName,'Eq_T4','dC','time',Double_FLAG=1)
    call NC_Def_Var(OutPath,CFM_OutName,'Eq_T34','dC','time',Double_FLAG=1)
    call NC_Def_Var(OutPath,CFM_OutName,'Eq_T3','dC','time',Double_FLAG=1)
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! Atm field of CF23
    IF(Clip_Run_OnOff == 0)THEN
    
    call NC_Create_File( OutPath,CFM_Atm_OutName, &
                         CFM_Par.ALon, CFM_Par.ALat,out_time_tick, &
                         'lon', 'lat','time', time_unit, 'X', Leap_Calendar_OnOff)
    
    call NC_Def_Var(OutPath,CFM_Atm_OutName,'KAtm','none','lon','time',Double_FLAG=1)
    call NC_Def_Var(OutPath,CFM_Atm_OutName,'RAtm','none','lon','time',Double_FLAG=1)
    
    call NC_Def_Var(OutPath,CFM_Atm_OutName,'IntDec','none','time',Double_FLAG=1);
    call NC_Def_Var(OutPath,CFM_Atm_OutName,'ItrAmp','none','time',Double_FLAG=1);
    
    call NC_Def_Var(OutPath,CFM_Atm_OutName,'Glob_Atm_U','m/s','lon','time',Double_FLAG=1)
    call NC_Def_Var(OutPath,CFM_Atm_OutName,'Glob_WWB_U','m/s','lon','time',Double_FLAG=1)
    
    ENDIF !!! IF(Clip_Run_OnOff == 0)
     
    write(*,*),'Chen-Fang 2023 Model OutPut Files Creating Completed'
    write(*,*),'------------------------------------------------------------------------------'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

END SUBROUTINE CFModel_Create_OutPut
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CFModel_Create_MonMe_OutPut(OutPath,CFM_Par,CFM_Stat,TimTick);

    character(len=*) :: OutPath
    TYPE(CFModel_Param) :: CFM_Par
    TYPE(CFModel_Stat)  :: CFM_Stat
    TYPE(TimeTick) :: TimTick  
    
    character(len=300) :: CFM_OutName
    character(len=300) :: CFM_Atm_OutName
    
    integer*8 :: IT
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
     real*8,allocatable :: out_time_tick(:) !!! tmp timetick size == OutPut_N_Event
     character(len=100) :: time_unit   !!! time_unit
     
     allocate(out_time_tick(MonMe_OutPut_N_Event)) !!! MonMe_OutPut_N_Event == MonMe_OutPut_N_Span
     
     !!! value for output_time_tick
     !!! because we have entering new month for TimTick
     DO  IT = 1, MonMe_OutPut_N_Event
         out_time_tick (IT) =  IT + (TimTick.Run_AbsIMon-1); 
         
         out_time_tick (IT) = out_time_tick(IT) * (365.0d0/12.0d0) 
     END DO
     
     !!! to prevent mixjudge of CDO (sub half-month)
     out_time_tick = out_time_tick - 0.5d0 * (365.0d0/12.0d0)
     
     !!! time_unit
     time_unit = &
     'days since '//trim(Math_YMD_To_Date(Start_CYear,Start_CMonth,Start_CDay))//' 00:00:0.0';

     !!! CFM_OutName
     CFM_OutName = 'MonMe_CF23_'//trim(TimTick.MonMe_Out_Date)//'_'//&
                   trim(Num2Str(MonMe_OutPut_Mon_Span))//'mon_R'//trim(Num2Str(TimTick.Run_Cycle_ID))//'.nc'

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     write(*,*),'----------------------------------------------------------------------------'
     write(*,*),'Creating Chen-Fang 2023 Model Monthly Mean OutPut Files ...'
     
     call NC_Create_File( OutPath,CFM_OutName, &
                          CFM_Par.OLon, CFM_Par.OLat,out_time_tick, &
                          'lon', 'lat','time', time_unit, 'X', Leap_Calendar_OnOff)
     
     call NC_Def_Var(OutPath,CFM_OutName,'KOcn','none','lon','time',Double_FLAG=1)
     call NC_Def_Var(OutPath,CFM_OutName,'ROcn','none','lon','time',Double_FLAG=1)
     call NC_Def_Var(OutPath,CFM_OutName,'Nd_Ocn_T','none','lon','time',Double_FLAG=1)
     
     call NC_Def_Var(OutPath,CFM_OutName,'Ocn_U','m/s','lon','time',Double_FLAG=1)
     call NC_Def_Var(OutPath,CFM_OutName,'Ocn_H','m','lon','time',Double_FLAG=1)
     call NC_Def_Var(OutPath,CFM_OutName,'Ocn_T','dC','lon','time',Double_FLAG=1)
     call NC_Def_Var(OutPath,CFM_OutName,'Atm_U','m/s','lon','time',Double_FLAG=1)
     
     
     call NC_Def_Var(OutPath,CFM_OutName,'ZAFK','dC/s','lon','time',Double_FLAG=1)
     call NC_Def_Var(OutPath,CFM_OutName,'THFK','dC/s','lon','time',Double_FLAG=1)
     call NC_Def_Var(OutPath,CFM_OutName,'RESQ','dC/s','lon','time',Double_FLAG=1)
     
     call NC_Def_Var(OutPath,CFM_OutName,'Eq_T4','dC','time',Double_FLAG=1)
     call NC_Def_Var(OutPath,CFM_OutName,'Eq_T34','dC','time',Double_FLAG=1)
     call NC_Def_Var(OutPath,CFM_OutName,'Eq_T3','dC','time',Double_FLAG=1)
     
     
     write(*,*),'Chen-Fang 2023 Model Monthly Mean OutPut Files Creating Completed'
     write(*,*),'------------------------------------------------------------------------------'
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END SUBROUTINE CFModel_Create_MonMe_OutPut
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CFModel_Write_OutPut(OutPath,CFM_Par,CFM_Stat,TimTick,WriteIT);

    character(len=*) :: OutPath
    TYPE(CFModel_Param) :: CFM_Par
    TYPE(CFModel_Stat)  :: CFM_Stat
    TYPE(TimeTick) :: TimTick  
    
    integer*4 :: WriteIT !!!! controlled by Ctrl_OutPut_Stream
    
    character(len=300) :: CFM_OutName
    character(len=300) :: CFM_Atm_OutName
    
    !!! CFM_OutName
    CFM_OutName = 'CF23_'//trim(TimTick.Out_Date)//'_'//trim(Num2Str(OutPut_Day_Span))//'day_'// &
                   trim(Num2Str(PutData_N_Day))//'day_R'//trim(Num2Str(TimTick.Run_Cycle_ID))//'.nc'
    
    CFM_Atm_OutName = 'CF23_Atm_'//trim(TimTick.Out_Date)//'_'//trim(Num2Str(OutPut_Day_Span))//'day_'// &
                       trim(Num2Str(PutData_N_Day))//'day_R'//trim(Num2Str(TimTick.Run_Cycle_ID))//'.nc'
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! BE CAREFUL : TimeTick.Run_AbsIT >> Res_AbsIT
    call NC_Write_One_IT(OutPath,CFM_OutName,'Res_AbsIT',dble(TimTick.Run_AbsIT),WriteIT);
    
    !!! CDay, CMonth, CDay into OutPut File
    call NC_Write_One_IT(OutPath,CFM_OutName,'CYear',dble(TimTick.CYear),WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'CMonth',dble(TimTick.CMonth),WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'CDay',dble(TimTick.CDay),WriteIT);
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    call NC_Write_One_IT(OutPath,CFM_OutName,'KOcn',CFM_Stat.KOcn,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'ROcn',CFM_Stat.ROcn,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'Nd_Ocn_T',CFM_Stat.Nd_Ocn_T,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'IntDec',CFM_Stat.IntDec,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'ItrAmp',CFM_Stat.ItrAmp,WriteIT);
    
    call NC_Write_One_IT(OutPath,CFM_OutName,'Ocn_U',CFM_Stat.Ocn_U,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'Ocn_H',CFM_Stat.Ocn_H,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'Ocn_T',CFM_Stat.Ocn_T,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'Atm_U',CFM_Stat.Atm_U,WriteIT);
    
    call NC_Write_One_IT(OutPath,CFM_OutName,'ZAFK',CFM_Stat.ZAFK,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'THFK',CFM_Stat.THFK,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'RESQ',CFM_Stat.RESQ,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'Ex_Ocn_Str',CFM_Stat.Ex_Ocn_Str,WriteIT);
    
    IF( CFModel_OutPut_2D_OnOff == 1) THEN
    
        call NC_Write_One_IT(OutPath,CFM_OutName,'Ocn_UU',CFM_Stat.Ocn_UU,WriteIT);
        call NC_Write_One_IT(OutPath,CFM_OutName,'Ocn_HH',CFM_Stat.Ocn_HH,WriteIT);
        call NC_Write_One_IT(OutPath,CFM_OutName,'Ocn_TT',CFM_Stat.Ocn_TT,WriteIT);
        call NC_Write_One_IT(OutPath,CFM_OutName,'Atm_UU',CFM_Stat.Atm_UU,WriteIT);
    
    END IF !!!! CFModel_OutPut_2D_OnOff
    
    
    call NC_Write_One_IT(OutPath,CFM_OutName,'Eq_T4',CFM_Stat.Eq_T4,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'Eq_T34',CFM_Stat.Eq_T34,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'Eq_T3',CFM_Stat.Eq_T3,WriteIT);
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    IF(Clip_Run_OnOff == 0)THEN
    !!!! Atm Field of CF23
    call NC_Write_One_IT(OutPath,CFM_Atm_OutName,'KAtm',CFM_Stat.KAtm,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_Atm_OutName,'RAtm',CFM_Stat.RAtm,WriteIT);
    
    call NC_Write_One_IT(OutPath,CFM_Atm_OutName,'IntDec',CFM_Stat.IntDec,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_Atm_OutName,'ItrAmp',CFM_Stat.ItrAmp,WriteIT);
    
    call NC_Write_One_IT(OutPath,CFM_Atm_OutName,'Glob_Atm_U',CFM_Stat.Glob_Atm_U,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_Atm_OutName,'Glob_WWB_U',CFM_Stat.Glob_WWB_U,WriteIT);
    
    ENDIF !!! IF(Clip_Run_OnOff == 0)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( Write_OutPut_Check == 1 ) THEN
    
    write(*,*),'------------------------------------------------------------------------------'
    write(*,*),'Writing Chen-Fang 2023 Model OutPut Files ...'
    
    write(6,101),WriteIT
    101 format(' WriteIT = ',I8);  
    write(*,*),'------------------------------------------------------------------------------'
    !!! Make sure not empty file output & check if openmp works correctly
    write(*,*),'Eq_T4       = '//trim(Num2Str(CFM_Stat.Eq_T4))//&
                '  Eq_T34       = '//trim(Num2Str(CFM_Stat.Eq_T34))//&
                '  Eq_T3        = '//trim(Num2Str(CFM_Stat.Eq_T3))//&
                '  Ex_Ocn_Str   = '//trim(Num2Str(CFM_Stat.Ex_Ocn_Str(1)));
    write(*,*),''
    write(*,*),'Chen-Fang 2023 Model OutPut Files Writing Completed'
    write(*,*),'------------------------------------------------------------------------------'
    
    END IF !!! Write_OutPut_Check 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
END SUBROUTINE CFModel_Write_OutPut
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CFModel_Write_MonMe_OutPut(OutPath,CFM_Par,CFM_Stat,TimTick,WriteIT);

    character(len=*) :: OutPath
    TYPE(CFModel_Param) :: CFM_Par
    TYPE(CFModel_Stat)  :: CFM_Stat
    TYPE(TimeTick) :: TimTick  
    
    integer*4 :: WriteIT !!!! controlled by Ctrl_OutPut_Stream
    
    character(len=300) :: CFM_OutName
    character(len=300) :: CFM_Atm_OutName
    
    !!! CFM_OutName
     CFM_OutName = 'MonMe_CF23_'//trim(TimTick.MonMe_Out_Date)//'_'//&
                   trim(Num2Str(MonMe_OutPut_Mon_Span))//'mon_R'//trim(Num2Str(TimTick.Run_Cycle_ID))//'.nc'
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    call NC_Write_One_IT(OutPath,CFM_OutName,'KOcn',CFM_Stat.MonMe_KOcn,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'ROcn',CFM_Stat.MonMe_ROcn,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'Nd_Ocn_T',CFM_Stat.MonMe_Nd_Ocn_T,WriteIT);
    
    call NC_Write_One_IT(OutPath,CFM_OutName,'Ocn_U',CFM_Stat.MonMe_Ocn_U,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'Ocn_H',CFM_Stat.MonMe_Ocn_H,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'Ocn_T',CFM_Stat.MonMe_Ocn_T,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'Atm_U',CFM_Stat.MonMe_Atm_U,WriteIT);
    
    call NC_Write_One_IT(OutPath,CFM_OutName,'ZAFK',CFM_Stat.MonMe_ZAFK,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'THFK',CFM_Stat.MonMe_THFK,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'RESQ',CFM_Stat.MonMe_RESQ,WriteIT);
    
    call NC_Write_One_IT(OutPath,CFM_OutName,'Eq_T4',CFM_Stat.MonMe_Eq_T4,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'Eq_T34',CFM_Stat.MonMe_Eq_T34,WriteIT);
    call NC_Write_One_IT(OutPath,CFM_OutName,'Eq_T3',CFM_Stat.MonMe_Eq_T3,WriteIT);
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( Write_OutPut_Check == 1 ) THEN
    
    write(*,*),'------------------------------------------------------------------------------'
    write(*,*),'Writing Chen-Fang 2023 Model Monthly Mean OutPut Files ...'
    
    write(6,101),WriteIT
    101 format(' WriteIT = ',I8);  
    write(*,*),'------------------------------------------------------------------------------'
    !!! Make sure not empty file output & check if openmp works correctly
    write(*,*),'Eq_T4       = '//trim(Num2Str(CFM_Stat.Eq_T4))//&
                '  Eq_T34       = '//trim(Num2Str(CFM_Stat.Eq_T34))//&
                '  Eq_T3        = '//trim(Num2Str(CFM_Stat.Eq_T3));
    write(*,*),''
    write(*,*),'Chen-Fang 2023 Model Monthly Mean OutPut Files Writing Completed'
    write(*,*),'------------------------------------------------------------------------------'
    
    END IF !!! Write_OutPut_Check 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
END SUBROUTINE CFModel_Write_MonMe_OutPut
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CFModel_Ctrl_ReStart_Stream (TimTick,CFM_Par,CFM_Stat) 

    TYPE(TimeTick)      :: TimTick
    TYPE(CFModel_Param) :: CFM_Par
    TYPE(CFModel_Stat)  :: CFM_Stat
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
    IF ( (TimTick.Res_FLAG_Count >= 1).AND.(mod(TimTick.Res_FLAG_Count, ReStart_N_Unit) == 0) )then
        
        !!! The following operations can be completed by ELBOM
        IF( CFModel_RunAlone_OnOff == 1) THEN 
            TimTick.Res_FLAG_Count = 0 !!!! very important, be careful
        ENDIF !!! CFModel_RunAlone_OnOff
    
        !!! update ReStart
        call CFModel_Write_ReStart(Restart_Path,CFM_Par,CFM_Stat,TimTick);
          
        IF (TimTick.RES_FLAG_Total_Count == ReStart_N_Unit * ReStart_Cycles) then
        
            IF ( TimTick.Run_AbsIT < TimTick.Run_AbsNT ) then 
            !!! to avoid STOP caused by ReStart at the end of one run
        
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!! The following operations can be completed by ELBOM
            IF( CFModel_RunAlone_OnOff == 1 ) then
            
            write(6,200) TimTick.Cyear,TimTick.Cmonth,TimTick.CDay, TimTick.Run_AbsIDay
            200 format(" Cyear = ",I4,"  Cmonth = ",I2,"  CDay = ",I2," Run_AbsIDay =",I10)
            write(*,*),''
            write(*,*),'ReStart Cycles Reached, Program Stopped'
            write(*,*),'------------------------------------------------------------------------------'
            write(6,100) Start_Cyear,Start_Cmonth,Start_CDay
            100  format(" Start_Cyear = ",I10,"  Start_Cmonth = ",I10,"  Start_CDay  = ",I10)
            write(6,101) End_Cyear,End_Cmonth,End_CDay
            101  format(" End_Cyear   = ",I10,"  End_Cmonth   = ",I10,"  End_CDay    = ",I10)
            write(6,102) Day_NT, TimTick.Run_AbsNT,TimTick.Run_AbsNDay
            102  format(" Day_NT      = ",I10,"  Run_AbsNT    = ",I10,"  Run_AbsNDay = ",I10)
            write(*,*)
            write(6,103) Leap_Calendar_OnOff
            103  format(' Leap_Calendar_OnOff   = ', I6);
            write(6,104) ReRun_OnOff
            104  format(' ReRun_Calendar_OnOff  = ', I6);
            write(6,105) ReStart_Cycles
            105  format(' ReStart_Cycles        = ', I6);
            write(6,106) ReStart_N_Unit
            106  format(' ReStart_N_Unit        = ', I6);
            write(*,*),'ReStart_DT_Unit       =  '//trim(ReStart_DT_Unit)
            write(*,*),'Show_TimeTick_DT_Unit =  '//trim(Show_TimeTick_DT_Unit)
            write(*,*),'Forcing_DT_Unit       =  '//trim(Forcing_DT_Unit)
            
            write(*,*),'-------------------------------------------------------------------------------'
            
            STOP !!! BE CAREFUL
            
            ENDIF !!! CFModel_RunAlone_OnOff 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            END IF !!! TimTick.Run_AbsIT < TimTick.Run_AbsNT
            
        END IF !!! TimTick.RES_FLAG_Total_Count == ReStart_N_Unit * ReStart_Cycles
    
    END IF !!! (TimTick.Res_FLAG_Count >= 1).AND.(mod(TimTick.Res_FLAG_Count, ReStart_N_Unit) == 0)

END SUBROUTINE CFModel_Ctrl_ReStart_Stream
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CFModel_Ctrl_OutPut_Stream(TimTick,CFM_Par,CFM_Stat) 

    TYPE(TimeTick)      :: TimTick
    TYPE(CFModel_Param) :: CFM_Par
    TYPE(CFModel_Stat)  :: CFM_Stat
    
    integer*4 :: WriteIT !!! which event to write (integer*4,BE CAREFUL)
    integer*8 :: AbsIT_resid !!! residual AbsIT
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( Inst_OutPut_OnOff == 1 ) THEN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! if reaches the first integration step of a new round of OutPut_Day_Span
    !!! create OutPut single file 
    IF ( Clip_Run_OnOff == 0 ) THEN
        IF ( mod (TimTick.Run_AbsIT, (Day_NT)*(OutPut_Day_Span) ) == 1 )THEN
    
            !!! update TimTick.Out_Date, help target the file to write
            TimTick.Out_Date = TimTick.Date
    
            call CFModel_Create_OutPut(OutPut_Path,CFM_Par,CFM_Stat,TimTick);
        
        ENDIF
    ELSEIF( Clip_Run_OnOff == 1 ) THEN
        IF ( TimTick.Run_AbsIT == (TimTick.Res_AbsIT + 1) ) THEN
        
            !!! update TimTick.Out_Date, help target the file to write
            TimTick.Out_Date = TimTick.Date
            
            call CFModel_Create_OutPut(OutPut_Path,CFM_Par,CFM_Stat,TimTick);
    
        END IF
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! if reaches the last integration step of a new round of PutData_N_Day
    !!! write OutPut single file
    IF ( Clip_Run_OnOff == 0 ) THEN
        IF ( (TimTick.Run_AbsIT>=1) .AND. (mod(TimTick.Run_AbsIT, (Day_NT)*(PutData_N_Day)) == 0) )THEN
    
             !!! be careful 
             AbsIT_resid = mod( TimTick.Run_AbsIT - 1, (Day_NT)*OutPut_Day_Span) + 1;
         
             !!! get writeIT
             WriteIT = int4( AbsIT_resid/((Day_NT)*(PutData_N_Day)) );
    
             call CFModel_Write_OutPut(OutPut_Path,CFM_Par,CFM_Stat,TimTick,WriteIT);
    
        END IF
    ELSEIF( Clip_Run_OnOff == 1 ) THEN
    
        IF ( (TimTick.Run_AbsIT>=1) .AND. (mod(TimTick.Run_AbsIT - TimTick.Res_AbsIT, (Day_NT)*(PutData_N_Day)) == 0) )THEN
        
            !!! be careful 
            WriteIT = int4 ( ceiling(dble(TimTick.Run_AbsIT - TimTick.Res_AbsIT)/dble(Day_NT*PutData_N_Day)));
        
            call CFModel_Write_OutPut(OutPut_Path,CFM_Par,CFM_Stat,TimTick,WriteIT);
    
        END IF
    
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END IF !!! IF ( Inst_OutPut_OnOff == 1 ) THEN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! Monthly Mean OutPut
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( MonMe_OutPut_OnOff == 1 ) THEN 
    
        !!! Create 
        IF( (TimTick.Arrive_CMonth_FLAG == 1).AND.( mod(TimTick.Run_AbsIMon, MonMe_OutPut_Mon_Span)==1 ) )THEN
    
            TimTick.MonMe_Out_Date = TimTick.Date
            
            call CFModel_Create_MonMe_OutPut(OutPut_Path,CFM_Par,CFM_Stat,TimTick);
        
        END IF
        
        !!! Write
        IF( TimTick.Reach_CMonth_FLAG == 1 )THEN
        
            WriteIT = int4( mod(TimTick.Run_AbsIMon-1, MonMe_OutPut_Mon_Span) +1 );
            
        
            call CFModel_Write_MonMe_OutPut(OutPut_Path,CFM_Par,CFM_Stat,TimTick,WriteIT);
        
        
        END IF 
    
    END IF !!! MonMe_OutPut_OnOff
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
END SUBROUTINE CFModel_Ctrl_OutPut_Stream
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CFModel_Recover_Field (TimTick,CFM_Par,CFM_Stat) 

    TYPE(TimeTick)      :: TimTick
    TYPE(CFModel_Param) :: CFM_Par
    TYPE(CFModel_Stat)  :: CFM_Stat
    
    integer*8 :: NOcn,NAtm,AI_NY,IX,IY
    
    real*8,allocatable :: tmp_Ka(:,:), tmp_Ra(:,:), tmp_Nd_Ocn_T(:,:)
    
    NOcn = CFM_Par.NOcn;
    NAtm = CFM_Par.NAtm;
    AI_NY = CFM_Par.AI_NY;
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(tmp_Ka(NAtm,1)); tmp_Ka = 0.0d0;
    allocate(tmp_Ra(NAtm,1)); tmp_Ra = 0.0d0;
    
    allocate(tmp_Nd_Ocn_T(NOcn,1)); 
    tmp_Nd_Ocn_T(1:NOcn,1) = CFM_Stat.Nd_Ocn_T(1:NOcn);
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    CFM_Stat.KOcn = CFM_Stat.UVec(1:NOcn,1);
    CFM_Stat.ROcn = CFM_Stat.UVec(NOcn+1:2*NOcn,1);
    CFM_Stat.Nd_Ocn_T = CFM_Stat.UVec(2*NOcn+1:3*NOcn,1);
    
    !!! Reconstructions (MATLAB code)
    !!! Hov_Ko_a = soln_store(1:No,:);
    !!! Hov_Ro_a = soln_store(No+1:2*No,:);
    !!! Hov_T_a = Dim_T * psi_0_eq * soln_store(2*No+1:3*No,:);
    !!!
    !!! Hov_U_a = Dim_U * ((Hov_Ko_a - Hov_Ro_a) * psi_0_eq + Hov_Ro_a/sqrt(2) * psi_2_eq);
    !!! Hov_H_a = Dim_H * ((Hov_Ko_a + Hov_Ro_a) * psi_0_eq + Hov_Ro_a/sqrt(2) * psi_2_eq);
    
    CFM_Stat.Ocn_T = ( CFM_Par.Dim_Ocn_T ) * ( CFM_Par.Eq_OPHI_0 ) * CFM_Stat.Nd_Ocn_T; 
    
    CFM_Stat.Ocn_U = CFM_Par.Dim_Ocn_U *  &
                    ( (CFM_Par.Eq_OPHI_0) * (CFM_Stat.KOcn - CFM_Stat.ROcn) +  &
                      (CFM_Par.Eq_OPHI_2/sqrt(2.0d0)) * (CFM_Stat.ROcn) ); 

    CFM_Stat.Ocn_H = CFM_Par.Dim_Ocn_H *  &
                    ( (CFM_Par.Eq_OPHI_0) * (CFM_Stat.KOcn + CFM_Stat.ROcn) +  &
                      (CFM_Par.Eq_OPHI_2/sqrt(2.0d0)) * (CFM_Stat.ROcn) ); 
    
    !!! get Eq_T3, Eq_T34, Eq_T4
    CFM_Stat.Eq_T3   = sum(CFM_Stat.Ocn_T(CFM_Par.T3_WInd : CFM_Par.T3_EInd))/ &
                       dble( CFM_Par.T3_EInd - CFM_Par.T3_WInd + 1 );
    
    CFM_Stat.Eq_T34  = sum(CFM_Stat.Ocn_T(CFM_Par.T34_WInd : CFM_Par.T34_EInd))/ &
                       dble( CFM_Par.T34_EInd - CFM_Par.T34_WInd + 1 );

    CFM_Stat.Eq_T4   = sum(CFM_Stat.Ocn_T(CFM_Par.T4_WInd : CFM_Par.T4_EInd))/ &
                       dble( CFM_Par.T4_EInd - CFM_Par.T4_WInd + 1 );
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IX)
    DO IX = 1,NOcn
    
       CFM_Stat.Ocn_UU(IX,1:AI_NY) =  CFM_Par.Dim_Ocn_U *  &
                                 ( (CFM_Par.OPHI_0) * (CFM_Stat.KOcn(IX) - CFM_Stat.ROcn(IX)) +  &
                                 (CFM_Par.OPHI_2/sqrt(2.0d0)) * (CFM_Stat.ROcn(IX)) );
       
       CFM_Stat.Ocn_HH(IX,1:AI_NY) =  CFM_Par.Dim_Ocn_H *  &
                                 ( (CFM_Par.OPHI_0) * (CFM_Stat.KOcn(IX) + CFM_Stat.ROcn(IX)) +  &
                                 (CFM_Par.OPHI_2/sqrt(2.0d0)) * (CFM_Stat.ROcn(IX)) );
       
       CFM_Stat.Ocn_TT(IX,1:AI_NY) =  &
       ( CFM_Par.Dim_Ocn_T ) * ( CFM_Par.OPHI_0 ) * CFM_Stat.Nd_Ocn_T(IX);
    
    END DO
    !$omp end parallel do 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! MATLAB code 
    !!! tsp_Ka_a = Ka_temp * Hov_T_a;
    !!! tsp_Ra_a = Ra_temp * Hov_T_a;
    !!! Hov_Ka_a = transpose(tsp_Ka_a);
    !!! Hov_Ra_a = transpose(tsp_Ra_a);
    !!! Hov_u_a = Dim_u * ((Hov_Ka_a - Hov_Ra_a) * phi_0_eq + Hov_Ra_a/sqrt(2) * phi_2_eq);
    !!! Hov_u_a = transpose(Hov_u_a);
    !!! Hov_u_a = Hov_u_a(1:No,:);
    
    tmp_Ka = matmul ( CFM_Par.Ka_temp, tmp_Nd_Ocn_T);
    tmp_Ra = matmul ( CFM_Par.Ra_temp, tmp_Nd_Ocn_T);
    
    
    CFM_Stat.Atm_U = CFM_Par.Dim_Atm_U * &  
                     ( (CFM_Par.Eq_APHI_0) * (tmp_Ka(1:NOcn,1) - tmp_Ra(1:NOcn,1)) +  &
                       (CFM_Par.Eq_APHI_2/sqrt(2.0d0)) * (tmp_Ra(1:NOcn,1)) ); 
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IX)
    DO IX = 1,NOcn
    
       CFM_Stat.Atm_UU(IX,1:AI_NY) = CFM_Par.Dim_Atm_U * & 
                                    ( (CFM_Par.APHI_0) * (tmp_Ka(IX,1) - tmp_Ra(IX,1)) +  &
                                      (CFM_Par.APHI_2/sqrt(2.0d0)) * (tmp_Ra(IX,1)) ); 
    
    END DO
    !$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! recover atm field
    
    CFM_Stat.KAtm = tmp_Ka(1:NAtm,1);
    CFM_Stat.RAtm = tmp_Ra(1:NAtm,1);
    
    !!! interannual wind 
    CFM_Stat.Glob_Atm_U = CFM_Par.Dim_Atm_U * &  
                          ( (CFM_Par.Eq_APHI_0) * (tmp_Ka(1:NAtm,1) - tmp_Ra(1:NAtm,1)) +  &
                            (CFM_Par.Eq_APHI_2/sqrt(2.0d0)) * (tmp_Ra(1:NAtm,1)) ); 
    
    !!! intraseasonal wind 
    CFM_Stat.Glob_WWB_U = CFM_Par.Dim_Atm_U * &
                          (CFM_Par.Eq_APHI_0) * (CFM_Par.glob_sp) * (CFM_Stat.ItrAmp); 
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    deallocate( tmp_Ka );
    deallocate( tmp_Ra );
    deallocate( tmp_Nd_Ocn_T );

END SUBROUTINE CFModel_Recover_Field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE IntDec_OneStep_Forward(TimTick,CFM_Par,CFM_Stat) 
    
    TYPE(TimeTick)      :: TimTick
    TYPE(CFModel_Param) :: CFM_Par
    TYPE(CFModel_Stat)  :: CFM_Stat
    
    real*8 :: tmp_randn, sgm_x
    integer*8 :: temp, size_sgm
    
    !!! get random number
    tmp_randn = CFM_Par.RandN( CFM_Stat.CFM_Run_AbsIT )
    size_sgm = size(CFM_Par.IdSgm);
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! interdecadal variability SODE (MATLAB code)
    !!! I = zeros(1,LL);I(1)=.5;
    !!! for i = 2:LL
    !!!    temp = round(I(i-1)*100) + 1;
    !!!    if temp<0
    !!!        temp = 1;
    !!!    elseif temp > n
    !!!        temp = n;
    !!!    end
    !!!    sgm_x = sgm(temp);
    !!!    I(i) = I(i-1) + (-lambda * (I(i-1) - m) * dt) + sgm_x * randn *sqrt(dt);
    !!! end
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    temp = nint(CFM_Stat.IntDec * 100.0d0 ) + 1;
    
    IF( temp < 0 ) THEN
        temp = 1
    ELSEIF( temp > size_sgm ) THEN
        temp = size_sgm
    ELSE
    END IF 
    
    sgm_x = CFM_Par.IdSgm(temp);
    
    !!! Euler Maruyama Scheme
    !!! I(i) = I(i-1) + (-lambda * (I(i-1) - m) * dt) + sgm_x * randn *sqrt(dt);
    CFM_Stat.IntDec = CFM_Stat.IntDec + &
             (-1.0d0) * CFM_Par.IdLambda * ( CFM_Stat.IntDec - CFM_Par.IdM ) * (CFM_Par.dt ) + &
             sgm_x * tmp_randn * sqrt( CFM_Par.dt );

    
END SUBROUTINE IntDec_OneStep_Forward
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CFModel_Update_EvoMat (TimTick,CFM_Par,CFM_Stat) 
    
    TYPE(TimeTick)      :: TimTick
    TYPE(CFModel_Param) :: CFM_Par
    TYPE(CFModel_Stat)  :: CFM_Stat
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer *8 :: NOcn,IX,IY,IT
    real*8 :: dt
    real*8 :: tmp_T4, tmp_phase,tmp_b
    
    NOcn   = CFM_Par.NOcn
    dt     = CFM_Par.dt
    tmp_phase  = (dt * CFM_Stat.CFM_Run_AbsIT * (34.0d0/365.0d0) + CFM_Par.phase )
    
    
    !!! update interannual + interdecadal evolution matrix
    
    !!! MATLAB code
    !Tw = mean(u0(west_node));
    !Te = mean(u0(east_node));
    !T4 = mean(u0(T4_node));
    CFM_Stat.T4 = sum(  CFM_Stat.UVec(2*NOcn + CFM_Par.T4_WInd : 2*NOcn + CFM_Par.T4_EInd, 1))/&
                  dble( CFM_Par.T4_EInd - CFM_Par.T4_WInd + 1);

    tmp_T4 = CFM_Stat.T4
    
    
    !% the M matrix depends on the decadal variability
    !for ss = 1:No
    !    M_31(ss,ss) = c1*(eta(ss)+I(i)*eta2(ss)*1.);
    !    M_32(ss,ss) = c1*(eta(ss)-I(i)*eta2(ss)*1.);
    !end
    !phase = -1/12; % seasonal cycle
    
    DO IX = 1,NOcn
    
       !!! M_31 diagonal element 
       CFM_Stat.EvoMat(2*NOcn + IX, IX ) = &
       (CFM_Par.c1)  * ( CFM_Par.eta(IX) + CFM_Stat.IntDec * CFM_Par.eta2(IX) * (1.0d0) );
       
       !!! M_32 digonal element
       CFM_Stat.EvoMat(2*NOcn + IX, NOcn + IX ) = &
       (CFM_Par.c1)  * ( CFM_Par.eta(IX) - CFM_Stat.IntDec * CFM_Par.eta2(IX) * (1.0d0) );
    
    END DO
    
    
    !% approximation of a cubic damping together with the seasonal effect
    !
    !b = c1 * zeta * alpha_q *(1.8- eta2/3+(0.2 + abs(mean(u0(T4_node))+0.4).* eta2 ).^2/5)...
    !    .* ( 1.+0.5*sin(2*pi*(dt*i*34/365+phase)) + 0.4*sin(2*pi*(dt*i*34/365+phase-1/12)) * eta2/4 - 0.15*sin(4*pi*(dt*i*34/365+phase-2/12))*eta/2.4);
    !
    !M_33 = -b .* eye(No) ;
    !
    !M = [M_11,M_12,M_13;
    !    M_21,M_22,M_23;
    !    M_31,M_32,M_33];
    
    DO IX = 1,NOcn
    
    !!! M_33 diagnoal element (BE VERY CAREFUL !!!)
    
    !!! c1 * zeta * alpha_q *(1.8- eta2/3+(0.2 + abs( mean(u0(T4_node)) + 0.4).* eta2 ).^2/5)...
    tmp_b = ( CFM_Par.c1 * CFM_Par.zeta * CFM_Par.alpha_q ) * & 
            ( 1.8d0 - CFM_Par.eta2(IX)/(3.0d0) + &
              (0.2d0) * ( 0.2d0 +  abs( tmp_T4 + 0.4d0 ) * CFM_Par.eta2(IX) )**(2.0d0)  );
    
    
    tmp_b = tmp_b * &
    ( 1.0d0 + &
      (         0.5d0          ) * ( 1.0d0 ) * sin ( 2.0d0 * pi * tmp_phase                      ) + &
      ( CFM_Par.eta2(IX)/4.0d0 ) * ( 0.4d0 ) * sin ( 2.0d0 * pi * ( tmp_phase - (1.0d0/12.0d0) ) ) + &
      ( CFM_Par.eta(IX)/2.4d0  ) * (-0.15d0) * sin ( 4.0d0 * pi * ( tmp_phase - (2.0d0/12.0d0) ) ) );
    
    
    CFM_Stat.EvoMat( 2*NOcn + IX, 2*NOcn + IX ) = (-1.0d0) * (tmp_b);
    

    END DO
    
    !!! Must Test the EvoMat value at first time step
    !!!write(*,*),'CFM_Run_AbsIT = ',CFM_Stat.CFM_Run_AbsIT
    !!!write(*,*),'tmp_T4        = ',tmp_T4
    !!!write(*,*),'tmp_b         = ',tmp_b
    !!!write(*,*),'M_31(end,end) = ',CFM_Stat.EvoMat(3*NOcn,NOcn)
    !!!write(*,*),'M_32(end,end) = ',CFM_Stat.EvoMat(3*NOcn,2*NOcn)
    !!!write(*,*),'M_33(end,end) = ',CFM_Stat.EvoMat(3*NOcn,3*NOcn)
    
END SUBROUTINE CFModel_Update_EvoMat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CFModel_Update_FrcVec (TimTick,CFM_Par,CFM_Stat) 
    
    TYPE(TimeTick)      :: TimTick
    TYPE(CFModel_Param) :: CFM_Par
    TYPE(CFModel_Stat)  :: CFM_Stat
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer *8 :: NOcn,IX,IY,IT
    real*8 :: dt, dt_noise
    real*8 :: tmp_T4, tmp_phase
    real*8 :: sgm_p
    
    NOcn     = CFM_Par.NOcn
    dt       = CFM_Par.dt
    dt_noise = dt/(10.0d0) 
    
    tmp_phase  = (dt * CFM_Stat.CFM_Run_AbsIT * (34.0d0/365.0d0) + CFM_Par.phase )
    
    CFM_Stat.T4 = sum(  CFM_Stat.UVec(2*NOcn + CFM_Par.T4_WInd : 2*NOcn + CFM_Par.T4_EInd, 1))/&
                  dble( CFM_Par.T4_EInd - CFM_Par.T4_WInd + 1);

    tmp_T4 = CFM_Stat.T4
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! MATLAB CODE
    !!!% noise is dependent on season; noise is also dependent on Walker
    !!!% circulation
    !!!sgm_p = (1.*tanh(1.*T4)+1.) *(1+0.6*cos(2*pi*(dt*i*34/365+phase-1/12)))* (1-I(i)*0.75)*1.6;
    !!!ap_mean = 0;
    !!!dp = 1/(365/34/12);% damping time of the wind bursts, now set to be 1 month.
    sgm_p = (1.0d0 * tanh(1.0d0*tmp_T4) + 1.0d0) * &
            (1.0d0 + 0.6d0 * cos(2.0d0 * pi *(tmp_phase - (1.0d0/12.0d0) ) ) ) * &
            (1.0d0 - CFM_Stat.IntDec * 0.75d0) * 1.6d0;
    
    
    !!!% Stochastic parameterization of the wind burst; its parameters depend on the state
    !!!ap_mean = 0;
    !!!for tt = 1:10
    !!!    ap_temp = ap_temp - dp * (ap_temp - ap_mean ) * dt_noise + sqrt(dt_noise) * sgm_p * rd_save(ii);
    !!!    ii = ii + 1;
    !!!end
    !!!ap(i) = ap_temp;
    
    DO IT = 1,10
       
       CFM_Stat.ItrAmp = CFM_Stat.ItrAmp + (-1.0d0)* CFM_Par.ItrDp * CFM_Stat.ItrAmp * dt_noise + &
       sqrt(dt_noise) * sgm_p * CFM_Par.Rd_Save( 10 * (CFM_Stat.CFM_Run_AbsIT-1) + IT)
    
    END DO
    
    
    !!!c2 = 0.1;
    !!!% add a mean forcing to compensate the flux from the cubic damping that causes non-zero mean 
    !!!f = [gamma * chi2 * c1 * ap(i) * sp / 2; - gamma * chi2 * c1 * ap(i) * sp / 3;...
    !!!    c1*eta2'* chi2*c2];
    
    CFM_Stat.FrcVec(1:NOcn,1) = ( CFM_Par.gamma * CFM_Par.chi2 * CFM_Par.c1 * CFM_Stat.ItrAmp ) * &
                                ( CFM_Par.sp /(2.0d0) );
    
    CFM_Stat.FrcVec(NOcn+1:2*NOcn,1) = ( CFM_Par.gamma * CFM_Par.chi2 * CFM_Par.c1 * CFM_Stat.ItrAmp ) * &
                                       ( CFM_Par.sp /(-3.0d0) );
    
    CFM_Stat.FrcVec(2*NOcn+1:3*NOcn,1) = (CFM_Par.c1 * CFM_Par.chi2 * CFM_Par.c2) * &
                                         (CFM_Par.eta2)
    
    
    !!! MUST TEST FrcVec value at first time step 
    !write(*,*),'CFM_Run_AbsIT    = ',CFM_Stat.CFM_Run_AbsIT
    !write(*,*),'sgm_p            = ',sgm_p
    !write(*,*),'CFM_Stat.ItrAmp  = ',CFM_Stat.ItrAmp
    !write(*,*),'F1(NOcn/2)       = ',CFM_Stat.FrcVec(NOcn   - NOcn/2,1)
    !write(*,*),'F2(NOcn/2)       = ',CFM_Stat.FrcVec(2*NOcn - NOcn/2,1)
    !write(*,*),'F3(NOcn/2)       = ',CFM_Stat.FrcVec(3*NOcn - NOcn/2,1)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! JYW ++ 
    !!! CFM_Stat.FrcVec = (CFM_Par.FrcVec_Amp) * CFM_Stat.FrcVec
    !!! JYW ++ 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END SUBROUTINE CFModel_Update_FrcVec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CFModel_Update_Ex_FrcVec (TimTick,CFM_Par,CFM_Stat) 
    
    TYPE(TimeTick)      :: TimTick
    TYPE(CFModel_Param) :: CFM_Par
    TYPE(CFModel_Stat)  :: CFM_Stat
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer *8 :: NOcn,IX,IY,IT
    real*8 :: dt, dt_noise
    
    NOcn     = CFM_Par.NOcn
    dt       = CFM_Par.dt
    dt_noise = dt/(10.0d0)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    CFM_Stat.Ex_FrcVec(1:NOcn,1) = ( CFM_Stat.Ex_Ocn_Str * CFM_Par.Re2Nd_Str )/(2.0d0);
    
    CFM_Stat.Ex_FrcVec(NOcn+1:2*NOcn,1) = ( CFM_Stat.Ex_Ocn_Str * CFM_Par.Re2Nd_Str )/(-3.0d0);
    
    CFM_Stat.Ex_FrcVec(2*NOcn+1:3*NOcn,1) = 0.0d0;
    
    
END SUBROUTINE CFModel_Update_Ex_FrcVec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CFModel_Get_MonMe_Stat(TimTick,CFM_Par,CFM_Stat)

    TYPE(TimeTick)      :: TimTick
    TYPE(CFModel_Param) :: CFM_Par
    TYPE(CFModel_Stat)  :: CFM_Stat
    
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer*8 :: NMode
    integer*8 :: IM 
    
    integer*8 :: NX,NY,NZ,ULNZ
    integer*8 :: IX,IY,IZ
    
    integer*8 :: Arrive_CMonth_FLAG !!! Arrive_FLAG: arrive at the start of month
    integer*8 :: Reach_CMonth_FLAG  !!! Reach_FLAG:  arrive at the end of month
    
    integer*8 :: CYear, CMonth


    Arrive_CMonth_FLAG = TimTick.Arrive_CMonth_FLAG
    Reach_CMonth_FLAG  = TimTick.Reach_CMonth_FLAG
    
    CYear  = TimTick.CYear
    CMonth = TimTick.CMonth
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !$omp parallel sections
    
    !!! MonMe OcnStat
    !$omp section
    call Update_MonMe_Var (CFM_Stat.KOcn, CFM_Stat.MonMe_KOcn, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    call Update_MonMe_Var (CFM_Stat.ROcn, CFM_Stat.MonMe_ROcn, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    call Update_MonMe_Var (CFM_Stat.Nd_Ocn_T, CFM_Stat.MonMe_Nd_Ocn_T, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    
    !$omp section
    call Update_MonMe_Var (CFM_Stat.Ocn_T, CFM_Stat.MonMe_Ocn_T, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    call Update_MonMe_Var (CFM_Stat.Ocn_H, CFM_Stat.MonMe_Ocn_H, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    call Update_MonMe_Var (CFM_Stat.Ocn_U, CFM_Stat.MonMe_Ocn_U, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    call Update_MonMe_Var (CFM_Stat.Atm_U, CFM_Stat.MonMe_Atm_U, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    !$omp section
    call Update_MonMe_Var (CFM_Stat.ZAFK, CFM_Stat.MonMe_ZAFK, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    call Update_MonMe_Var (CFM_Stat.THFK, CFM_Stat.MonMe_THFK, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    call Update_MonMe_Var (CFM_Stat.RESQ, CFM_Stat.MonMe_RESQ, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    !$omp section
    call Update_MonMe_Var (CFM_Stat.Eq_T4,  CFM_Stat.MonMe_Eq_T4, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    call Update_MonMe_Var (CFM_Stat.Eq_T34, CFM_Stat.MonMe_Eq_T34, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    call Update_MonMe_Var (CFM_Stat.Eq_T3,  CFM_Stat.MonMe_Eq_T3, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    
    !$omp end parallel sections 
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE CFModel_Get_MonMe_Stat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CFModel_OneStep_Forward(TimTick,CFM_Par,CFM_Stat) 
    
    TYPE(TimeTick)      :: TimTick
    TYPE(CFModel_Param) :: CFM_Par
    TYPE(CFModel_Stat)  :: CFM_Stat
    
    real*8 :: Mon2dt, Sec2dt
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer *8 :: NOcn,IX,IY,IT
    real*8 :: dt
    
    NOcn   = CFM_Par.NOcn
    dt     = CFM_Par.dt
    Mon2dt = (365.0d0/12.0d0) * (2.0d0) * dt
    Sec2dt = ( Mon2dt/(365.0d0/12.0d0) )/(86400.0d0)
    
    !!! update evolution matrix
    call CFModel_Update_EvoMat (TimTick,CFM_Par,CFM_Stat) 
    
    !!! update westerly wind burst (FrcVec)
    call CFModel_Update_FrcVec (TimTick,CFM_Par,CFM_Stat)
    
    
    !!! update external forcing vec (Ex_FrcVec)
    call CFModel_Update_Ex_FrcVec (TimTick,CFM_Par,CFM_Stat)
    
    
    !!! record heat budget tendency (Step One)
    CFM_Stat.TEND =  ( CFM_Par.Dim_Ocn_T ) * ( CFM_Par.Eq_OPHI_0 ) * CFM_Stat.Nd_Ocn_T;
    
    CFM_Stat.ZAFK =  ( CFM_Par.Eq_OPHI_0 * CFM_Par.Dim_Ocn_T * CFM_Par.c1 ) * (CFM_Stat.IntDec) * &
                       CFM_Par.eta2 * (CFM_Stat.KOcn - CFM_Stat.ROcn); 
    CFM_Stat.ZAFK = CFM_Stat.ZAFK * Sec2dt
    
    CFM_Stat.THFK =  ( CFM_Par.Eq_OPHI_0 * CFM_Par.Dim_Ocn_T * CFM_Par.c1 ) * &
                       CFM_Par.eta * (CFM_Stat.KOcn + CFM_Stat.ROcn); 
    CFM_Stat.THFK = CFM_Stat.THFK * Sec2dt
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! main solver 
    !!! u = u0 + M*u0*dt  + f * dt;
    !!! u0 = u;
    CFM_Stat.UVec = CFM_Stat.UVec + dt * matmul(CFM_Stat.EvoMat, CFM_Stat.UVec) 
    
    
    CFM_Stat.UVec = CFM_Stat.UVec + (CFM_Par.FrcVec_Amp) * dt * CFM_Stat.FrcVec
    
    
    CFM_Stat.UVec = CFM_Stat.UVec + (CFM_Par.Ex_FrcVec_Amp) * dt * CFM_Stat.Ex_FrcVec;
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!! record heat budget tendency (step two)
    CFM_Stat.TEND = &
    (  ( CFM_Par.Dim_Ocn_T ) * ( CFM_Par.Eq_OPHI_0 ) * CFM_Stat.UVec(2*NOcn+1:3*NOcn,1) - CFM_Stat.TEND )/dt;
    
    !!! dC/dt* ==> dC/s
    CFM_Stat.TEND = CFM_Stat.TEND * Sec2dt
    
    CFM_Stat.RESQ = CFM_Stat.TEND - CFM_Stat.ZAFK - CFM_Stat.THFK;
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! interdecadal must be behind 
    !!! interdecadal variability SODE
    call IntDec_OneStep_Forward (TimTick,CFM_Par,CFM_Stat) 
    
    !!! Reconstruct (Recover) Field 
    call CFModel_Recover_Field (TimTick,CFM_Par,CFM_Stat)
    
    !!! MUST TEST UVec value at first time step
    !!! write(*,*),'CFM_Run_AbsIT    = ',CFM_Stat.CFM_Run_AbsIT
    !!! write(*,*),'UVec1(NOcn/2)    = ',CFM_Stat.UVec(NOcn   - NOcn/2,1)
    !!! write(*,*),'UVec2(NOcn/2)    = ',CFM_Stat.UVec(2*NOcn - NOcn/2,1)
    !!! write(*,*),'UVec3(NOcn/2)    = ',CFM_Stat.UVec(3*NOcn - NOcn/2,1)
    
    
END SUBROUTINE CFModel_OneStep_Forward
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CFModel_Get_Eq_Nudging_Val (GrdInfo,TimTick,CFM_Par,CFM_Stat,Clim_OcnStat,OcnStat) 
   
    !!! give equatorial nudging to the model of ELBOM 
    TYPE(GridInfo)  :: GrdInfo
    TYPE(TimeTick)      :: TimTick
    TYPE(CFModel_Param) :: CFM_Par
    TYPE(CFModel_Stat)  :: CFM_Stat
    TYPE(OceanStat) :: OcnStat
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer*8 :: IX,IY,IT
    integer*8 :: NX,NY
    real*8 :: N4_Skew_A, N4_Skew_B, N4_Skew_C, N4_Skew_M
    real*8 :: N3_Skew_A, N3_Skew_B, N3_Skew_C, N3_Skew_M
    real*8 :: Skew_SSTA_N4, Skew_SSTA_N3
    
    !!! VERY IMPORTANT for ReStart
    OcnStat.Eq_Nudging_T4 = CFM_Stat.Eq_T4;
    OcnStat.Eq_Nudging_T34 = CFM_Stat.Eq_T34;
    OcnStat.Eq_Nudging_T3 = CFM_Stat.Eq_T3;
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    N4_Skew_A = Clim_OcnStat.N4_Skew_A
    N4_Skew_B = Clim_OcnStat.N4_Skew_B
    N4_Skew_C = Clim_OcnStat.N4_Skew_C
    N4_Skew_M = Clim_OcnStat.N4_Skew_M
    
    N3_Skew_A = Clim_OcnStat.N3_Skew_A
    N3_Skew_B = Clim_OcnStat.N3_Skew_B
    N3_Skew_C = Clim_OcnStat.N3_Skew_C
    N3_Skew_M = Clim_OcnStat.N3_Skew_M
    
    !!! Skew_SSTA_N4 = psgn_skew(CFM_Stat.Eq_T4,N4_Skew_A,N4_Skew_B,N4_Skew_C,N4_Skew_M)
    !!! Skew_SSTA_N3 = psgn_skew(CFM_Stat.Eq_T3,N3_Skew_A,N3_Skew_B,N3_Skew_C,N3_Skew_M)

    Skew_SSTA_N4 = CFM_Stat.Eq_T4
    Skew_SSTA_N3 = CFM_Stat.Eq_T3
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! Get (Skewed) Regression SST
    OcnStat.Reg_SST = (Clim_OcnStat.SSTA_Reg_N4) * (Skew_SSTA_N4) + &
                      (Clim_OcnStat.SSTA_Reg_N3) * (Skew_SSTA_N3) + &
                      (Clim_OcnStat.SSTA_Reg_Const)
    
    !!! Fill Galapagos islands (288x64 grid)
    !!! OcnStat.Reg_SST(261:262,32) = (0.5d0) * ( OcnStat.Reg_SST(261:262,31) + OcnStat.Reg_SST(261:262,33) )
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!! linear interpolation of Ocn_T from CFM_stat (natural boundary condition)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! !$omp parallel do schedule (static) &
    !!! !$omp default(shared) &
    !!! !$omp private(IX)
    !!! DO IX = 1,(GrdInfo.NX + 2)
    !!!
    !!!    OcnStat.Eq_Nudging_SST(IX) = Math_FindValue_In_MZSeq &
    !!!    ( MZseq = CFM_Par.OLon, Fseq = CFM_Stat.Ocn_T, Coord = GrdInfo.PLon(IX), Tone = 'A')
    !!!
    !!! OcnStat.Eq_Nudging_SST(IX) = Nudging_SSTA_Amp * OcnStat.Eq_Nudging_SST(IX)
    !!!
    !!! END DO !!! IX
    !!! !$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IX,IY)
    DO IX = 1,(GrdInfo.NX + 2)
        DO IY = 1,(GrdInfo.NY + 2)
    
        !!! OcnStat.Nudging_SST(IX, IY) =  Nudging_SSTA_Amp * OcnStat.Eq_Nudging_SST(IX) * CFM_Par.OPHI_0(IY)
        
        !!! OcnStat.Nudging_SST(IX,IY) = OcnStat.Eq_Nudging_SST(IX)
        
        OcnStat.Nudging_SST(IX,IY) = Nudging_SSTA_Amp * OcnStat.Reg_SST(IX,IY)

        !!! OcnStat.Nudging_SST(IX,IY) = psgn(OcnStat.Nudging_SST(IX,IY))

    
        END DO
    END DO
    !$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! write(*,*),'OcnStat.PLon (NX/2)           = ',GrdInfo.PLon(GrdInfo.NX/2) 
    !!! write(*,*),'OcnStat.Eq_Nudging_SST (NX/2) = ',OcnStat.Eq_Nudging_SST(GrdInfo.NX/2) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
END SUBROUTINE CFModel_Get_Eq_Nudging_Val
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
END MODULE MO_CF23MODEL
