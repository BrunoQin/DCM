MODULE MO_TYPEDEF
    
USE MO_NCTOOLS
USE MO_MATHTOOLS
USE ISO_C_BINDING
USE torch_wrapper

IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! include 'ocean_comblk.f90'
!!!
!!! InPut/ReStart/OutPut Path
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!!! character(len = 100), parameter :: Model_Path = 'D:/MyModels/My_ELBOM/CF23_ICM/CF23_ICM_Fortran/'
character(len = 100), parameter :: Model_Path = '../'
character(len = 100), parameter :: Tsub_model_loc = "/home/jywang/FTA/lib/offline_cesm2_online_soda23_tsub.pt"
character(len = 100), parameter :: Tau_model_loc = "/home/jywang/FTA/lib/offline_reg_obs_tau.pt"

character(len=100),parameter :: ReStart_Path = trim(Model_Path)//'restart/';
character(len=100),parameter :: OutPut_Path = trim(Model_Path)//'output/';

character(len=100),parameter :: Grid_Path = trim(Model_Path)//'input/';
character(len=100),parameter :: ClimOcn_Path = trim(Model_Path)//'input/Model_Clim_OceanStat/';
character(len=100),parameter :: Forcing_Path = trim(Model_Path)//'input/Model_Forcing/';

character(len=100),parameter :: EigMode_FName = 'WOA18_EigMode.nc';
character(len=100),parameter :: Clim_PsTauX_FName = 'UGrid_MonClim_PsTauX_1991_2020.nc';
character(len=100),parameter :: Clim_PsTauY_FName = 'VGrid_MonClim_PsTauY_1991_2020.nc';
character(len=100),parameter :: Clim_TauX_FName   = 'UGrid_MonClim_TauX.nc';
character(len=100),parameter :: Clim_TauY_FName   = 'VGrid_MonClim_TauY.nc'
character(len=100),parameter :: CMIP_Clim_TauX_FName   = 'CMIP_UGrid_MonClim_TauX.nc';
character(len=100),parameter :: CMIP_Clim_TauY_FName   = 'CMIP_VGrid_MonClim_TauY.nc';

character(len=100),parameter :: PsTauX_Forcing_FName = 'UGrid_PsTauX_Forcing.nc';
character(len=100),parameter :: PsTauY_Forcing_FName = 'VGrid_PsTauY_Forcing.nc';
character(len=100),parameter :: TauX_Forcing_FName   = 'UGrid_TauX_Forcing.nc';
character(len=100),parameter :: TauY_Forcing_FName   = 'VGrid_TauY_Forcing.nc';
character(len=100),parameter :: CMIP_TauX_Forcing_FName   = 'LP_CMIP_UGrid_TauX_Forcing.nc';
character(len=100),parameter :: CMIP_TauY_Forcing_FName   = 'LP_CMIP_VGrid_TauY_Forcing.nc';

character(len=100),parameter :: OBS_Nudging_FName = 'GODAS_Nudging_N3N4.nc';

character(len=100),parameter :: Clim_OceanStat_FName = 'PUVGrid_MonClim_OceanStat.nc';
character(len=100),parameter :: EOF_Linear_Tau_Model_FName = 'EOF_Linear_Tau_Model.nc';
character(len=100),parameter :: Reg_Linear_Tau_Model_FName = 'Reg_Linear_Tau_Model.nc';
character(len=100),parameter :: CFModel_Param_FName = 'CF23_Param_1000yr.nc';


!!! test TSUB-AI
character(len=100),parameter :: SSH_Test_FName='MonAn_Mask_SSH_CESM2_piControl_0001_0010.nc';
character(len=100),parameter :: TSUB_Test_FName='MonAn_Mask_TEMP_50m_CESM2_piControl_0001_0010.nc';

!!! CNOP MJO Perturbation
character(len=100),parameter :: MJO_Forcing_EEOF_FName='CF23_ENBOM_MJO_PsTauXY_Forcing_EEOF.nc'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! InPut/ReStart/OutPut Stream Namelist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! OpenMP Speed-UP
integer*4,parameter :: MY_OMP_THREADS = 24

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Two ways to run : (Default) Calendar_Run 
!!! Calendar_Run : start & end timetick assigned by Calendar
!!! Clip_Run     : read intial value from ReStart File, Run_NDay (not Abs but Relative) 
!!!                directly assigned             

!!! start timetick calendar year,month,day for one run 
!!! start from 00:00:00
integer*8,parameter :: Start_CYear = 1, Start_CMonth = 1, Start_CDay = 1 
!!! end at 24:00:00
integer*8,parameter :: End_CYear   = 20, End_CMonth = 12, End_CDay = 31

!!! calendar 365day (0) or not (1)
integer*8,parameter  :: Leap_Calendar_OnOff = 0; 

!!! time interval to update RESTART file: DT_Name = 'D' (1 day),'M' (month), 'Y' (year)
!!! if reach the end integration step of calendar year/month/day
character,parameter :: ReStart_DT_Unit = 'Y'
integer*8,parameter :: ReStart_N_Unit = 10

!!! STOP after number restart_cycles
integer*8,parameter :: ReStart_Cycles = 500; !!! default == 1

!!! show timetick calender when reaching one day ('D') or  one month ('M') or one year ('Y')
character,parameter :: Show_TimeTick_DT_Unit = 'M' 

!!! show writeIT and time series value 
integer*8,parameter :: Write_OutPut_Check = 0

!!!! SIMUTAENOUS OUTPUT 
!!!! create single OutPut file every (OutPut_Day_Span) day 
!!!                and just at the very first integration step of a new round (Day_Span)
integer*8,parameter :: OutPut_Day_Span = 365*20

!!!! putdata into OutPut file every (PutData_N_Day) day
!!!!               and just at the end integration step of a new round (N_Day)
integer*8,parameter :: PutData_N_Day   = 365

!!!! how many events each OutPut file has, mod (OutPut_Day_Span, PutData_N_Day) should == 0 !!!
integer*8,parameter :: OutPut_N_Event  = ceiling(dble(OutPut_Day_Span)/dble(PutData_N_Day))


!!!! create single MonMe OutPut file every (OutPut_Mon_Span) months
!!!                and just at the very first integration step of a new round (Mon_Span)
integer*8,parameter :: MonMe_OutPut_Mon_Span = 12*20;

!!!! putdata into MonMe OutPut file every (PutData_N_Mon) Month
!!!!               and just at the end integration step of a new round (N_Month)
integer*8,parameter :: MonMe_PutData_N_Mon   = 1;

!!!! how many events each MonMe OutPut file has
integer*8,parameter :: MonMe_OutPut_N_Event  =  MonMe_OutPut_Mon_Span;

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! the time resolution of external forcing data 
!!!! monthly forcing ('M') 
!!!! daily forcing ('D') 
!!!! clim  forcing ('C') 
!!!! free coupled run ('F')
!!!! direct nudged assimilation ('A')
!!!! indirect wind nudged assimilation ('N')
character :: Forcing_DT_Unit = 'F'

!!! nudging onoff
integer*8 :: SSTA_Nudging_OnOff = 1;          !!! nudging of SSTA
integer*8 :: Wind_Nudging_OnOff = 1;          !!! nudging of wind stresses  
integer*8 :: Add_MJO_Forcing_OnOff = 0;       !!! Add MJO_Forcing_OnOff

real*8    :: Add_MJO_EastWind_Ratio = 0.0d0;  !!! Add MJO EastWind Ratio (0.0d0~1.0d0) 
integer*8 :: Add_ENBOM_MJO_Forcing_OnOff = 1; !!! Add MJO_Forcing onto ENBOM
integer*8 :: Add_CF23_MJO_Forcing_OnOff  = 1; !!! Add MJO_Forcing onto CF23 

integer*8 :: Eq_Direct_Assign   = 0;     !!! direct assign value when nudging

real*8 :: SSTA_Nudging_Scale = 60.0d0    !!! unit: day
real*8 :: Wind_Nudging_Ratio = 1.0d0     !!! the ratio of external wind

real*8 :: SSTA_Nudging_Radius = 30.0d0   !!! unit: degree

!!! SSTA nudging
real*8,parameter :: Nudging_SSTA_Amp  = 1.0d0;

!!! SSTA nudging range (degree) 
real*8,parameter :: Nudging_WLon = 30.0d0; !!! 117.0d0;
real*8,parameter :: Nudging_ELon = 292.0d0;
real*8,parameter :: Nudging_SLat = -30.0d0;
real*8,parameter :: Nudging_NLat = +30.0d0;

!!! west wind bias (N/m^2) 
real*8,parameter :: MJO_WestWind_Bias = 0.0000d0
real*8,parameter :: PsTauXY_To_KelCoeff = 1.8828d0 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! in MAIN PROGRAM  
integer*8,parameter :: Run_Cycles    = 1;  !!! maximum cycles of running the main program (not main cycle) 
integer*8,parameter :: Run_OnOff     = 1;   !!! run main cycle (1) or not (0)
integer*8,parameter :: ReRun_OnOff   = 0;   !!! restart run (1) or not (0)
integer*8,parameter :: Solver_OnOff  = 1;   !!! use OceanSolver (1) or not (0)
integer*8,parameter :: OutPut_OnOff  = 1;   !!! create OutPut stream (1) or not (0)
integer*8,parameter :: ReStart_OnOff = 1;   !!! create ReStart stream (1) or not (0)

!!! INSTANT OUTPUT 
integer*8,parameter :: Inst_OutPut_OnOff = 0; 

!!! MONTHLY MEAN OUTPUT (FOR LONG-TERM SIMULATION TO SAVE SPACE)
integer*8,parameter :: MonMe_OutPut_OnOff = 1; 


!!! Clip_Run on (1) or not (0, Calendar Run)
integer*8,parameter :: Clip_Run_OnOff = 0
integer*8,parameter :: Clip_Run_NDay = 365*2

integer*8,parameter:: CFModel_RunAlone_OnOff = 0;   !!! CF23 Model Run Alone (1) or not (0)
integer*8,parameter :: CFModel_OutPut_2D_OnOff = 0; !!! CF23 Model OutPut 2D Field (1) or not (0) to save space

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! solve SSTA field (1) or not (0)
integer*8,parameter :: Solve_SSTA_OnOff = 1;

!!! solve TSUB field (1) or not (0) using AI
integer*8,parameter :: Solve_TSUB_AI_OnOff = 1;

!!! solve TAUXY field by TAUXY_AI (1), TAUXY_AI_EOF(2), TAUXY_EOF(3), TAUXY_REG(0)
integer*8,parameter :: Solve_TAUXY_AI_OnOff = 0;

!!! output SSTA budget (1) or not (0)
integer*8,parameter :: OutPut_HeatBDG_OnOff  = 0;

!!! output simple SSTA budget (1) or not (0) (ZAFK + THFK + RESFK + NUDGE)
integer*8,parameter :: OutPut_Simple_HeatBDG_OnOff = 0;

!!! output upper layer 4D field (1) or not (0)
integer*8,parameter :: OutPut_UpperLayer_OnOff = 0; 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! time constant
real*8,parameter    :: Ocn_NumDT   = 3600.0d0 * 3.0d0;  
real*8,parameter    :: CFModel_NumDT   = 3600.0d0 * 12.0d0;
integer*8,parameter :: CFModel_SplitNT = ceiling( CFModel_NumDT / Ocn_NumDT ) 

!!! Tsub model coupled N Day
integer*8,parameter :: Tsub_Couple_N_Day = 5
integer*8,parameter :: Tau_Couple_N_Day  = 5

!!! ocean model integration DT (second: unit)
real*8,parameter    :: Day_DT = 86400.0d0; !!! second each day
real*8,parameter    :: Year_DT = Day_DT * (365.0d0); !!! second each year
integer*8,parameter :: Day_NT = ceiling(Day_DT/Ocn_NumDT);
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! Physical Parameter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Truncation mode number of linear EOF Tau Model
integer*8,parameter :: Tau_EOF_Trunc_NMode = 24

!!! Trunction mode number of tsub recon
integer*8,parameter :: Tsub_EOF_Recon_OnOff = 1
integer*8,parameter :: Tsub_EOF_Trunc_NMode = 96

!!!! Correction Mask of Tsub
integer*8,parameter :: Tsub_Correct_OnOff = 0
!!!! real*8,parameter :: MaxVal_TSUB_OutPut_Mask = 1.5d0  !!! maxvalue of modulation mask of TSUB

!!! Tsub Mix up with SST
integer*8,parameter :: Tsub_Slope_Mix_OnOff = 0      !!! Mix up with SST 
real*8,parameter :: Tsub_Slope_Mix_Amp      = 1.0d0  !!! Minmum Mix Ratio of TSUB with SST

      
real*8,parameter :: Wstr_Proj_Amp    = 1.9d0    !!! wind stress projection amplitude factor 
real*8,parameter :: Couple_Strength  = 1.0d0    !!! (1.0d0) !!! air-sea coupled strength
real*8,parameter :: Couple_TAUX_Bias = 0.0d0    !!! (-0.00062d0)  !!! couple taux bias (N/m^2)

real*8,parameter :: TSUB_InPut_Amp   = 1.0d0    !!! TSUB input (SSH) amplitude factor
real*8,parameter :: TSUB_OutPut_Amp  = 0.9d0    !!! TSUB output amplitude factor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real*8,parameter :: EFold_C0  = 60.0d0;  !!! E-folding scale (m/s) for TAUXY_AI in/output

real*8,parameter :: PhaseLocking_OnOff = 0.0d0;   !!! (0.11d0) phase-locking sin fuction
real*8,parameter :: Phase_Difference   = 0.0d0; !!! (10.0d0) phase difference of sin function

real*8, parameter :: DynDamp_OnOff  = 1.0d0; !!! damping and dissipative effect for OceanSolver

!!! for OutPut Upper Layer 4D Field Depth cut off 
real*8,parameter :: Ref_UppLayer_Depth = 300.0d0; 
!!! reference upper layer depth (m), usually 300m or 150m

real*8,parameter :: RhoWat = 1000.0d0; !!!! reference sea water density (kg/m^3)

!!! bulk formula for wind stress
real*8,parameter :: CoeffDrag = 1.25d-3; !!! drag coefficient 
real*8,parameter :: RhoAtm = 1.23d0;     !!! air density (kg/m^3)

!!! horizontal momentum & mass viscosity coefficient 
real*8, parameter :: XY_Diff_Co_UVcur = (DynDamp_OnOff) * 3.0d3  !!! 3.0d3 (m^2/s)  
real*8, parameter :: XY_Diff_Co_Pres = (DynDamp_OnOff) * 1.0d3   !!! 1.0d3 (m^2/s)  

!!! vertical momentum & mass diffusion coefficient: the A of  (A/C^2)
real*8, parameter :: Z_Diff_Co_UVcur = (DynDamp_OnOff) * 5.0d-8 !!! (m^2/s^3)  
real*8, parameter :: Z_Diff_Co_Pres  = (DynDamp_OnOff) * 5.0d-8 !!! (m^2/s^3) 

!!! Ekman Layer damping 
real*8, parameter :: Ekm_UVShr_Damp = (1.0d0)/(2.0d0 * Day_DT); 

!!! control amplitude of nine feedbacks of SSTA
real*8,parameter :: UATBX_Amp = 1.0d0,  UBTAX_Amp = 1.0d0,  UATAX_Amp = 1.0d0
real*8,parameter :: VATBY_Amp = 1.0d0,  VBTAY_Amp = 1.0d0,  VATAY_Amp = 1.0d0
real*8,parameter :: WATBZ_Amp = 1.0d0,  WBTAZ_Amp = 0.75d0, WATAZ_Amp = 0.75d0 

!!! SSTA heat diffusion X & Y direction
real*8,parameter :: DiffX_SSTA = 3.0d3, DiffY_SSTA = 3.0d3, DiffZ_SSTA = 0.0d0 

!!! SSTA Damping
real*8,parameter :: Damp_SSTA = (1.0d0)/(100.0d0 * Day_DT);
real*8,parameter :: SSTA_DampPower = 1.0d0; !!! linear damping (1) or cubic damping (3)

!!! Linear Regression Bias for Nino4 & Nino3 (degree C)
real*8,parameter :: Reg_N4_Bias = 0.0d0;
real*8,parameter :: Reg_N3_Bias = 0.0d0;
real*8,parameter :: Reg_Couple_TAUX_Bias = 0.0d0; !!! (N/m^2) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TYPE BoundIndex
    
     integer*8 :: SInd, NInd, WInd, EInd
    
END TYPE BoundIndex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TYPE GridInfo

    integer*8 :: NX,NY,NZ,NMode,ClimNT
    real*8 :: ReDX,ReDY,ReDZ
    real*8 :: DLon,DLat
    
    real*8,allocatable :: Depth(:)
    real*8,allocatable :: IMode(:)
    
    real*8,allocatable :: ULDepth(:) !!!! upper layer depth array 
    integer*8 :: ULNZ !!! upper layer numbers
    
    real*8,allocatable :: PLon(:),PLat(:),ULon(:),ULat(:),VLon(:),VLat(:)
    real*8,allocatable :: PReX(:),PReY(:),UReX(:),UReY(:),VReX(:),VReY(:)
    real*8,allocatable :: PMask(:,:),UMask(:,:),VMask(:,:)
    
    !!! e-folding Pmask to preprocess input for TAUXY_AI
    real*8,allocatable :: Tau_PMask(:,:) 
    
    real*8 :: ReDT !!! ocean model time step
    
    !!! Nino Region Boundary Index
    TYPE(BoundIndex) :: Nino4_Bd, Nino3p4_Bd, Nino3_Bd
    integer*8 :: Nino4_PointNum, Nino3p4_PointNum, Nino3_PointNum
    
    !!! Equatorial Region Boundary Index
    TYPE(BoundIndex) :: Eq_Nino4_Bd, Eq_Nino3p4_Bd, Eq_Nino3_Bd
    integer*8 :: Eq_Nino4_PointNum, Eq_Nino3p4_PointNum, Eq_Nino3_PointNum
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Sponge Boundary Condition Alpha
    !!!
    !!! Un+1* = (1-alp)*Un+1 + alp*U0
    !!! (Un+1*-Un+1)/DT = alp(U0-Un+1)/DT 
    !!! (Un+1*-Un)/DT = (1-alp)*F(Un) + alp*(U0-Un)/DT
    !!!
    !!! alp = (1-(IY)/(Num))^(Pow)
    !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer*8 :: Sponge_Num = 6 !!! = 9 !!! how many sponge layers 
    integer*8 :: Sponge_Pow = 3 !!! power function
    
    real*8,allocatable :: Sponge_Alp(:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END TYPE GridInfo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

TYPE TimeTick

    !!! current time tick of the ocean state vector 
    
    integer*8 :: Run_Cycle_ID  !!! the ID of repeatedly running the whole program 

    integer*8 :: Res_AbsIT   !!! absolute integration IT for RESTART in one run
    
    integer*8 :: Run_AbsIT   !!! absolute integration IT for RUN in one run 
    integer*8 :: Run_AbsIDay !!! absolute integration IDay for RUN in one run 
    integer*8 :: Run_AbsIMon !!! absolute integration IMon for RUN in one run 
    
    integer*8 :: Run_AbsNT   !!! absolute integration NT in one run
    integer*8 :: Run_AbsNDay !!! absolute integration NDay in one run
    
    !!! Calendar year, month, day
    integer*8 :: CYear,CMonth,CDay
    
    !!! day sum in a year (1 to 365 or 366)
    integer*8 :: YrDaySum
    
    !!! Reach FLAG & Arrive_FLAG
    !!! if Reach_FLAG ==1, reached the end of day/month/year
    integer*8 :: Reach_CDay_FLAG, Reach_CMonth_FLAG, Reach_CYear_FLAG
    
    !!! if Arrive_FLAG == 1, arrived at the start of day/month/year 
    integer*8 :: Arrive_CDay_FLAG, Arrive_CMonth_FLAG, Arrive_CYear_FLAG
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! Restart FLAG Count 
    integer*8 :: Res_FLAG_Count       !!! change from 0 to ReStart_N_Unit
    
    !!! Restart FLAG Count (only increase)
    integer*8 :: Res_FLAG_Total_Count !!! if Total_Count == ReStart_Cycles, PROGRAM STOP
    
    !!! Date
    !!! CDO standard unit: 'days since 1981-01-01 00:00:0.0'
    character(len=100) :: Date   !!! Date for TimeTick: YYYY-MM-DD 
    
    !!! Helping to target Write single OutPut file Name (Date)
    character(len=100) :: Out_Date 
    !!! (Carve on gunwale of a moving boat)
    
    !!! Helping to target Write single Monthly Mean OutPut 
    character(len=100) :: MonMe_Out_Date
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    TYPE(ftorchmodel) :: Tsub_model  !!! tsub pytorch model pointer
    TYPE(ftorchmodel) :: Tau_model   !!! tau  pytorch model pointer
    
    INTEGER(C_INT) :: use_gpu = 1  !!! should be compatible with your script module device
    
    
END TYPE TimeTick

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

TYPE OceanStat
     
     !!!! vertical baroclinic mode coffecient for Ucur,Vcur,Pres
     real*8,allocatable :: Co_Ucur(:,:,:),Co_Vcur(:,:,:),Co_Pres(:,:,:)
     
     !!!! real physical field for upper-layer (300m) Ucur,Vcur,DzDt
     real*8,allocatable :: UL_Ucur(:,:,:),UL_Vcur(:,:,:),UL_DzDt(:,:,:)
     
     !!!! Ekman Layer Shear
     real*8,allocatable :: Ekm_Ushr(:,:),Ekm_Vshr(:,:)
     
     !!! surface field associated with SSTA tendency equation
     real*8,allocatable :: Usrf(:,:),Vsrf(:,:),Wsrf(:,:),SSH(:,:),SST(:,:),TSUB(:,:)
     
     !!! Usrf, SSH, TSUB, SST averaged over Nino4, Nino3
     real*8 :: SSHA_N4, SSHA_N34, SSHA_N3
     real*8 :: UA_N4,   UA_N34,   UA_N3
     
     real*8 :: TSUB_N4, TSUB_N34, TSUB_N3
     real*8 :: SSTA_N4, SSTA_N34, SSTA_N3
     real*8 :: Eq_SSTA_N4, Eq_SSTA_N34, Eq_SSTA_N3
     
     !!! couple with Chen-Fang 2023 Model (multimodel communication)
     real*8,allocatable :: Eq_Nudging_SST(:)
     real*8,allocatable :: Nudging_SST(:,:)
     real*8 :: Eq_Nudging_T4, Eq_Nudging_T34, Eq_Nudging_T3
     
     !!! regression model output
     real*8,allocatable :: Reg_SST(:,:)
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     !!! monthly mean state of anomalies for output analysis
     real*8,allocatable :: MonMe_Usrf(:,:), MonMe_Vsrf(:,:), MonMe_Wsrf(:,:)
     real*8,allocatable :: MonMe_SSH(:,:), MonMe_SST(:,:),  MonMe_TSUB(:,:)
     
     real*8 :: MonMe_SSTA_N4, MonMe_SSTA_N34, MonMe_SSTA_N3
     real*8 :: MonMe_Eq_Nudging_T4, MonMe_Eq_Nudging_T34, MonMe_Eq_Nudging_T3


END TYPE OceanStat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

TYPE Clim_OceanStat

     !!! eigen mode (NZ,NMODE)
     real*8,allocatable ::  EigMode(:,:,:) !!! Baroclinic Eigen Mode
     real*8,allocatable ::  EigC(:,:)      !!! Equivalent Kevlin Speed  
     real*8,allocatable ::  EigC_Pow2(:,:) !!! The (EigC)^2
     
     real*8,allocatable ::  Co_Wstr(:,:)   !!! coefficient for wind stress curve 
     !!! Co_Wstr * (WindStress)/ (RefWat) = wind force on single mode
     
     !!! IZ_EigMode & DZ_EigMode are on PUV grid, not DzDt grid 
     !!! IZ_EigMode for recovering DzDt on PUV grid
     real*8,allocatable :: IZ_EigMode(:,:,:) !!! Vertical Integral of Eigen Mode
     !!! DZ_EigMode for recovering Rho on PUV grid
     real*8,allocatable :: DZ_EigMode(:,:,:) !!! Vertical Derivative of Eigen Mode
     
     !!! clim pseudo-stress
     real*8,allocatable :: PsTauXBar(:,:,:) 
     real*8,allocatable :: PsTauYBar(:,:,:) 
     
     !!! clim wind-stress
     real*8,allocatable :: TauXBar(:,:,:) 
     real*8,allocatable :: TauYBar(:,:,:) 
     
     !!! CMIP clim wind stress
     real*8,allocatable :: CMIP_TauXBar(:,:,:)
     real*8,allocatable :: CMIP_TauYBar(:,:,:)
     
     !!! clim Usrf,Vsrf,Wsrf,SSH,SST,MLD (mixed layer depth)
     real*8,allocatable :: UsrfBar(:,:,:),VsrfBar(:,:,:) !!! UGrid,VGrid
     real*8,allocatable :: WsrfBar(:,:,:),SSHBar(:,:,:),SSTBar(:,:,:),TSUBBar(:,:,:) !!! PGrid
     real*8,allocatable :: TCDBar(:,:,:), MLDBar(:,:,:) !!! PGrid
     
     !!! CMIP Clim SSH, SST 
     real*8,allocatable :: CMIP_SSHBar(:,:,:)
     real*8,allocatable :: CMIP_SSTBar(:,:,:)
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     !!! Traditional Linear EOF PC Learning  Tau Model
     integer*8 :: Tau_EOF_NMode
     real*8,allocatable :: SSH_SST_scale(:) 
     real*8,allocatable :: SSH_SST_EOF(:,:,:,:)   !!! two level (NX+2,NY+2,2,EOF_NMode)
     real*8,allocatable :: TAUX_TAUY_EOF(:,:,:,:) !!! two level (NX+2,NY+2,2,EOF_NMode)
     
     !!! EOF linear regression coefficient
     real*8,allocatable :: Tau_AMat_12Mon(:,:,:) !!! (EOF_NMode,EOF_NMode,12)
     real*8,allocatable :: Tau_BVec_12Mon(:,:)   !!! (EOF_NMode,12)
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     !!! Nino Regression (REG) coefficient matrix (:,:,12Mon)
     real*8,allocatable :: TAUX_Reg_N4(:,:,:),TAUX_Reg_N3(:,:,:),TAUX_Reg_Const(:,:,:)
     real*8,allocatable :: TAUY_Reg_N4(:,:,:),TAUY_Reg_N3(:,:,:),TAUY_Reg_Const(:,:,:)
     
     !!! 2D regession
     real*8,allocatable :: SSTA_Reg_N4(:,:), SSTA_Reg_N3(:,:),SSTA_Reg_Const(:,:)
     
     !!! skewed 2D Regression
     real*8,allocatable :: Skew_SSTA_Reg_N4(:,:),Skew_SSTA_Reg_N3(:,:),Skew_SSTA_Reg_Const(:,:)
     
     !!! skewed tanh function param
     real*8 :: N4_Skew_A, N4_Skew_B, N4_Skew_C, N4_Skew_M 
     real*8 :: N3_Skew_A, N3_Skew_B, N3_Skew_C, N3_Skew_M 

     !!! standard deviation
     real*8,allocatable :: Std_SSTA(:,:)

     !!! Thermocline FeedBack Mask
     real*8,allocatable :: THFK_Mask(:,:)
     
     !!! TSUB OutPut Mask (12 months)
     real*8,allocatable :: TSUB_OutPut_Mask(:,:,:)
     
     !!! EOF TSUB to remove noise of TSUB-AI
     real*8,allocatable :: EOF_TSUB(:,:,:)
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     !!!! CNOP MJO Perturbation
     !!!! NEEOF = 40, NYear = 40 (1981-2020), NOcn(CF23Model)=56
     integer*8 :: NEEOF = 40
     integer*8 :: CFM_NOcn = 56
     
     real*8,allocatable :: EEOF_PsTauX(:,:,:,:) !!! (NX+2,NY+2,365,NEEOF)
     real*8,allocatable :: EEOF_PsTauY(:,:,:,:) !!! (NX+2,NY+2,365,NEEOF)
     real*8,allocatable :: PC_PsTauXY(:,:)      !!! (NEEOF,NYear)
     real*8,allocatable :: EEOF_KelCoeff(:,:,:) !!! (56,365,NEEOF)
     
     
END TYPE Clim_OceanStat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

TYPE AllFrc_OceanStat !!! right-hand forcing terms F(U)

    real*8,allocatable :: AllFrc_Co_Pres(:,:,:)
    real*8,allocatable :: AllFrc_Co_Ucur(:,:,:)
    real*8,allocatable :: AllFrc_Co_Vcur(:,:,:)
    
    !!! SSTA tendency equation feedbacks
    !!! Advective feedback: UATBX, Ekman feedback: WATBZ, Thermocline feedback : WBTAZ
    real*8,allocatable :: UATBX(:,:), UBTAX(:,:), UATAX(:,:)
    real*8,allocatable :: VATBY(:,:), VBTAY(:,:), VATAY(:,:)
    real*8,allocatable :: WATBZ(:,:), WBTAZ(:,:), WATAZ(:,:) 
    real*8,allocatable :: SSTA_RESQ(:,:) !!!! residual heat flux Q
    
    !!! simple HeatBDG feedbacks
    real*8,allocatable :: ZAFK(:,:)  !!! zonal advective feedback == UATBX
    real*8,allocatable :: THFK(:,:)  !!! thermocline feedback     == WBTAZ
    real*8,allocatable :: RESFK(:,:) !!! residual model feedback  == RESFK
    real*8,allocatable :: NUDGE(:,:) !!! external nudging term    == NUDGE
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! monthly mean state of anomalies
    real*8,allocatable :: MonMe_ZAFK(:,:),MonMe_THFK(:,:),MonMe_RESFK(:,:),MonMe_NUDGE(:,:)

END TYPE AllFrc_OceanStat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

TYPE SurfaceFlux !!! sea surface flux into the ocean 

    !!! momentum flux (wind stress)
    real*8,allocatable :: TAUX(:,:) !!! (NX+1,NY)
    real*8,allocatable :: TAUY(:,:) !!! (NX,NY+1)
    
    !!! for offline test of tau model
    real*8,allocatable :: TAUX_Test(:,:) !!! (NX+1,NY)
    real*8,allocatable :: TAUY_Test(:,:) !!! (NX,NY+1)
    
    real*8 :: TAUX_N4, TAUX_N34, TAUX_N3
    real*8 :: TAUX_N4_Test, TAUX_N34_Test, TAUX_N3_Test

    real*8,allocatable :: SSH(:,:)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! monthly mean state of anomalies
    real*8,allocatable :: MonMe_TAUX(:,:), MonMe_TAUY(:,:)
    real*8,allocatable :: MonMe_TAUX_Test(:,:), MonMe_TAUY_Test(:,:)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! CNOP MJO Perturbation
    integer*8 :: NEEOF = 40
    integer*8 :: CFM_NOcn = 56
    
    real*8,allocatable :: MJO_PsTauX(:,:)   !!! (NX+1,NY)
    real*8,allocatable :: MJO_PsTauY(:,:)   !!! (NX,NY+1)
    real*8,allocatable :: PC_MJO_PsTauXY(:) !!! (NEEOF)
    real*8,allocatable :: MJO_KelCoeff(:)   !!! (CFM_NOcn)
    

END TYPE SurfaceFlux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

TYPE SurfaceFlux_Hist !!! sea surface flux history (usually one year from external data)

    real*8,allocatable :: TAUX_Hist(:,:,:) !!! (NX+1,NY,:)
    real*8,allocatable :: TAUY_Hist(:,:,:) !!! (NX,NY+1,:)
    
    real*8,allocatable :: Nino4_Hist(:)
    real*8,allocatable :: Nino34_Hist(:)
    real*8,allocatable :: Nino3_Hist(:)

    real*8,allocatable :: SSH_Hist(:,:,:) 


END TYPE SurfaceFlux_Hist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TYPE CFModel_Param

    !!! non-dimensional time dt (12.00 hour) 
    !!! real*8 :: dt = (1.0d0/60.0d0) * ( (365.0d0)/(34.0d0 * 12.0d0) ); !!! calendar: 360_day
    real*8 :: dt = (360.0d0/365.0d0) * &
                   (1.0d0/60.0d0) * ( (365.0d0)/(34.0d0 * 12.0d0) ); !!! calendar: 365_day

    !!! non-dimensional dx
    real*8 :: dx = (8.0d0/3.0d0)/(128.0d0)
 
    real*8 :: LaDim = 40000.0d0; !!! global (km)
    real*8 :: LoDim = 17500.0d0; !!! Pacific (km)
    real*8 :: La = (8.0d0/3.0d0)
    real*8 :: Lo = (17500.d0/40000.0d0) * (8.0d0/3.0d0) 

    !!!! dimensional parameter
    integer*8 :: NAtm = 128 !!! number of amosphere
    integer*8 :: NOcn = 56  !!! number of ocean

    !!! ocean boundary
    !!! real*8 :: OWBD = 117.0d0, OEBD = 282.0d0;
    real*8 :: OWBD = 126.0d0, OEBD = 280.6875d0;
    real*8 :: OSBD = -31.5d0, ONBD = 31.5d0;

    integer*8 :: AI_NX = 288
    integer*8 :: AI_NY = 64

    !!! equivalent barolcinic ocean & atm speed (m/s)
    real*8 :: CAtm = 50.0d0
    real*8 :: COcn = 2.5d0;

    !!! dimensional factor
    real*8 :: Dim_Ocn_U = 0.5d0  !!! (m/s)
    real*8 :: Dim_Ocn_H = 20.8d0 !!! (m)
    real*8 :: Dim_Ocn_T = 1.5d0  !!! (dC)
    real*8 :: Dim_Atm_U = 5.0d0  !!! (m/s)

    real*8 :: chi1 = 0.65d0/(2.0d0); !!! projection from ocn to atm
    real*8 :: chi2 = 0.65d0*(2.0d0); !!! projection from atm to ocn

    !!! T3_node = No*2+31:No*2+53;
    !!! T4_node = No*2+14:No*2+31;
    !!! T34_node = No*2+25:No*2+43;
    integer*8 :: T3_WInd  = 31, T3_EInd  = 53
    integer*8 :: T4_WInd  = 14, T4_EInd  = 31
    integer*8 :: T34_WInd = 25, T34_EInd = 43

    !!! ocean,atm (lon,lat) 
    real*8,allocatable :: OLon(:),OLat(:),OReX(:),OReY(:)
    real*8,allocatable :: ALon(:),ALat(:),AReX(:),AReY(:)

    !!! non-dimensional OReY,AReY ==> ONdReY,ANdReY
    real*8,allocatable :: ONdReY(:),ANdReY(:)

    !!! meridional parabolic cylinder function mode (zero order + second order)
    real*8,allocatable :: OPHI_0(:),OPHI_2(:)
    real*8,allocatable :: APHI_0(:),APHI_2(:)

    !!! meridional mode's value on the equator
    real*8 :: Eq_OPHI_0, Eq_OPHI_2
    real*8 :: Eq_APHI_0, Eq_APHI_2

    !!! Ocn_T to Atm Kelvin & Rossby matrix
    real*8,allocatable :: Ka_temp(:,:),Ra_temp(:,:) !!! (128 x 56, NAtm x NOcn)

    !!! random number series (Time Step: 1/(365*2) Year = 12.0d0 hour )
    real*8,allocatable :: RandN(:)   !!! for interdecadal variability
    real*8,allocatable :: Rd_Save(:) !!! for westerly wind burst

    !!! Interdecadal variability (m, lamda, sgm(100)) in MATLAB
    real*8 :: IdM, IdLambda
    real*8,allocatable :: IdSgm(:) 

    !!! MAIN SOLVER (Interannual + Interdecadal)
    !!! take reduced gravity 0.03, mean thermocline depth 50, typical ocean velocity 0.5; Fr = U/sqrt(gh); 
    real*8 :: c1 = 0.15d0

    !!! Thermocline coefficient eta depends on x
    !!! Zonal advective coefficient eta2 depends on x
    real*8,allocatable :: eta(:),eta2(:) 
    real*8,allocatable :: Lin_M(:,:) !!! Linear Evolution Matrix (NOcn x 3, NOcn x 3)

    real*8 :: zeta = 8.7d0; !!! latent heating exchange coefficient
    real*8 :: phase = (-1.0d0)/12.0d0; !!! seasonal cycle

    !!! qc = 7; % latent heating multiplier coefficient
    !!! qe = 0.09296; % latent heating exponential coefficient
    !!! Tbar = 25/1.5; %1.5 is dimension, mean SST
    !!! tauq = 15; % latent heating adjustment rate
    !!! alpha_q = qc*qe*exp(qe*Tbar)/tauq; % latent heating factor
    real*8 :: alpha_q = ( (7.0d0) * (0.09296d0) * exp( (0.09296d0) * (25.0d0/1.5d0) ) )/(15.0d0)


    !!! Westerly Wind Burst noise
    !!! damping time of the wind bursts, now set to be 1 month
    real*8 :: ItrDp = (1.0d0)/( (365.0d0)/( 34.d0 * 12.0d0) ); 

    real*8,allocatable :: sp(:) !!! westerly wind burst profile
    real*8,allocatable :: glob_sp(:) !!!! global westerly burst profile

    real*8 :: gamma = 6.529d0   !!! wind stress coefficient

    !!! add a mean forcing to compensate the flux from the cubic damping that causes non-zero mean
    real*8 :: c2 = 0.1d0

    !!! FrcVec_Amp + Ex_FrcVec_Amp
    real*8 :: FrcVec_Amp = 1.0d0;
    real*8 :: Ex_FrcVec_Amp = 1.0d0;


    real*8 :: CFbeta = 2.28d-11;
    real*8 :: delta = 0.1d0; !!! long-wave scaling factor
    real*8 :: delta_Ocn = 0.1d0; !!! arbitary constant 

    real*8 :: reduced_g = 0.03d0;
    real*8 :: H_Ocn = (2.5d0*2.5d0)/(0.03d0); !!! (COcn*COcn)/reduced_g
    real*8 :: rho_Ocn = 1000.0d0;

    !!! very important
    !!! Dim_Str = delta * sqrt(CFbeta/CAtm) * (H_Ocn * rho_Ocn) * (COcn * COcn) * delta_Ocn;
    real*8 :: Dim_Str = &
              (0.1d0) * sqrt(2.28d-11/50.0d0) * (1000.0d0 * (2.5d0 * 2.5d0/0.03d0 ) ) * (2.5d0 * 2.5d0) * (0.1d0);

    !!! Re2NdStr = (c1) * (1./Dim_Str);
    real*8 :: Re2Nd_Str = (0.15d0)/ &
              ( (0.1d0) * sqrt(2.28d-11/50.0d0) * (1000.0d0 * (2.5d0 * 2.5d0/0.03d0 ) ) * (2.5d0 * 2.5d0) * (0.1d0) );
    

END TYPE CFModel_Param
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TYPE CFModel_Stat

    !!! 12.0d0 Hour IT
    integer*8 :: CFM_Run_AbsIT = 0 

    !!! Nondimenisonal Kelvin,Rossby Coefficient
    real*8,allocatable :: KAtm(:),RAtm(:) !!! (NAtm)
    real*8,allocatable :: KOcn(:),ROcn(:),Nd_Ocn_T(:) !!! (NOcn), Nd_Ocn_T: nom-dimensional Ocn_T 

    !!! state col vector (3*NOcn,1) [KOcn; ROcn; Ocn_T/(Dim_Ocn_T * OPHI_0_Eq)]
    !!! Nd_Ocn_T == Ocn_T/(Dim_Ocn_T * OPHI_0_Eq)
    real*8,allocatable :: UVec(:,:)  

    !!! 1D equatorial field (NOcn)
    real*8,allocatable :: Ocn_U(:) !!! ocean current anomaly
    real*8,allocatable :: Ocn_H(:) !!! ocean thermocline depth anomaly
    real*8,allocatable :: Ocn_T(:) !!! ocean temperature anomaly (SSTA)
    real*8,allocatable :: Atm_U(:) !!! interannual atmospheric wind anomaly 

    !!! 2D Field after reconstruction using parabolic cylinder function (NOcn,AI_NY)
    real*8,allocatable :: Ocn_UU(:,:)
    real*8,allocatable :: Ocn_HH(:,:) 
    real*8,allocatable :: Ocn_TT(:,:)            
    real*8,allocatable :: Atm_UU(:,:)

    !!! Nino T4,T34,T3 value on the euqator
    real*8 :: Eq_T4, Eq_T34, Eq_T3
    real*8 :: T4 !!! T4 = mean(u0(T4_node)); 

    !!! interdecadal variability (satisfy state-dependent stochastic ODE)
    real*8 :: IntDec

    !!! intraseasonal variability amplitude (ap in MATLAB)
    real*8 :: ItrAmp

    !!! state-dependent evolution Matrix ( Un+1 = M(Un) * Un ) 
    real*8,allocatable :: EvoMat(:,:) !!! (NOcn x3 , NOcn x 3)

    !!! state-dependent random forcing vec F(Un)
    real*8,allocatable :: FrcVec(:,:) !!! (NOcn x 3, 1)

    !!! (MJO) external forcing vec
    real*8,allocatable :: Ex_FrcVec(:,:) !!! (NOcn x 3,1)


    real*8,allocatable :: Glob_Atm_U(:) !!! global interannual atmospheric wind anomaly (NAtm)
    real*8,allocatable :: Glob_WWB_U(:) !!! global westerly wind burst wind anaomaly (NAtm)

    real*8,allocatable :: Ex_Ocn_Str(:) !!! External Zonal Ocean Wind Stress (dimensional N/m^2)
                                        !!! or projection coefficient of oceanic Kelvin wave

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!! simplified heat budget of CF23 (unit: dC/s not dC/dt*, be careful)
    real*8,allocatable :: TEND(:) !!! Total Tendency
    real*8,allocatable :: ZAFK(:) !!! Zonal Advective FeedBack
    real*8,allocatable :: THFK(:) !!! Thermocline FeedBack
    real*8,allocatable :: RESQ(:) !!! Residual Heating Term

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!! monthly mean state of anomalies
    real*8,allocatable :: MonMe_KOcn(:)     
    real*8,allocatable :: MonMe_ROcn(:)    
    real*8,allocatable :: MonMe_Nd_Ocn_T(:)

    real*8,allocatable :: MonMe_Ocn_U(:) 
    real*8,allocatable :: MonMe_Ocn_H(:) 
    real*8,allocatable :: MonMe_Ocn_T(:) 
    real*8,allocatable :: MonMe_Atm_U(:) 

    real*8,allocatable :: MonMe_TEND(:)
    real*8,allocatable :: MonMe_ZAFK(:)
    real*8,allocatable :: MonMe_THFK(:)
    real*8,allocatable :: MonMe_RESQ(:)
    
    real*8 :: MonMe_Eq_T4, MonMe_Eq_T34, MonMe_Eq_T3
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END TYPE CFModel_Stat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   
END MODULE MO_TYPEDEF
