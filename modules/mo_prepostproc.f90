MODULE MO_PREPOSTPROC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Pre-Processing & Post-Processing ToolBox
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE MO_NCTOOLS
USE MO_MATHTOOLS
USE MO_TYPEDEF
USE MO_OCEANSOLVER

IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Alloc_OceanStat(GrdInfo, OcnStat);

    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat)  :: OcnStat
    integer*8 :: NX,NY,NZ,NMode
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    NZ = GrdInfo.NZ
    NMode = GrdInfo.NMode
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! coefficient of baroclinic mode
    allocate(OcnStat.Co_Ucur(NX+1,NY,NMode));   OcnStat.Co_Ucur = 0.0d0;
    allocate(OcnStat.Co_Vcur(NX,NY+1,NMode));   OcnStat.Co_Vcur = 0.0d0;
    allocate(OcnStat.Co_Pres(NX+2,NY+2,NMode)); OcnStat.Co_Pres = 0.0d0;
    
    !!! upper field
    allocate(OcnStat.UL_Ucur(NX+1,NY,GrdInfo.ULNZ));      OcnStat.UL_Ucur = 0.0d0;
    allocate(OcnStat.UL_Vcur(NX,NY+1,GrdInfo.ULNZ));      OcnStat.UL_Vcur = 0.0d0;
    allocate(OcnStat.UL_DzDt(NX+2,NY+2,GrdInfo.ULNZ));    OcnStat.UL_DzDt = 0.0d0;
    
    !!! ekman layer shear
    allocate(OcnStat.Ekm_Ushr(NX+1,NY)); OcnStat.Ekm_Ushr = 0.0d0;
    allocate(OcnStat.Ekm_Vshr(NX,NY+1)); OcnStat.Ekm_Vshr = 0.0d0;
    
    !!! surface field
    allocate(OcnStat.Usrf(NX+1,NY));   OcnStat.Usrf = 0.0d0;
    allocate(OcnStat.Vsrf(NX,NY+1));   OcnStat.Vsrf = 0.0d0;
    allocate(OcnStat.Wsrf(NX+2,NY+2)); OcnStat.Wsrf = 0.0d0;
    allocate(OcnStat.SSH(NX+2,NY+2));  OcnStat.SSH  = 0.0d0;
    allocate(OcnStat.SST(NX+2,NY+2));  OcnStat.SST  = 0.0d0;
    allocate(OcnStat.TSUB(NX+2,NY+2)); OcnStat.TSUB = 0.0d0;
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    OcnStat.SSHA_N4 = 0.0d0;  OcnStat.SSHA_N34 = 0.0d0;  OcnStat.SSHA_N3 = 0.0d0;
    OcnStat.UA_N4   = 0.0d0;  OcnStat.UA_N34   = 0.0d0;  OcnStat.UA_N3   = 0.0d0;
    OcnStat.TSUB_N4 = 0.0d0;  OcnStat.TSUB_N34 = 0.0d0;  OcnStat.TSUB_N3 = 0.0d0;
    OcnStat.SSTA_N4 = 0.0d0;  OcnStat.SSTA_N34 = 0.0d0;  OcnStat.SSTA_N3 = 0.0d0;
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! couple with Chen-Fang 2023 Model (multimodel communication)
    allocate(OcnStat.Eq_Nudging_SST (NX+2) ); 
    OcnStat.Eq_Nudging_SST = 0.0d0;
    
    allocate(OcnStat.Nudging_SST (NX+2,NY+2) );
    OcnStat.Nudging_SST = 0.0d0;
    
    allocate(OcnStat.Reg_SST (NX+2,NY+2) );
    OcnStat.Reg_SST = 0.0d0;
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    allocate( OcnStat.MonMe_Usrf(NX+1,NY) ); OcnStat.MonMe_Usrf = 0.0d0;
    allocate( OcnStat.MonMe_Vsrf(NX,NY+1) ); OcnStat.MonMe_Vsrf = 0.0d0;
    
    allocate( OcnStat.MonMe_Wsrf(NX+2,NY+2) ); OcnStat.MonMe_Wsrf = 0.0d0;
    allocate( OcnStat.MonMe_SSH(NX+2,NY+2) );  OcnStat.MonMe_SSH  = 0.0d0;
    allocate( OcnStat.MonMe_SST(NX+2,NY+2) );  OcnStat.MonMe_SST  = 0.0d0;
    allocate( OcnStat.MonMe_TSUB(NX+2,NY+2) ); OcnStat.MonMe_TSUB = 0.0d0;
    
    
END SUBROUTINE Alloc_OceanStat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Alloc_TimeTick (TimTick)

    TYPE(TimeTick) :: TimTick
    CHARACTER(100),TARGET :: Tsub_modloc, Tau_modloc
    
    !!! intialization
    TimTick.Res_AbsIT = 0; 
    
    TimTick.Run_AbsIT = 0;
    TimTick.Run_AbsIDay = 0;
    
    TimTick.CYear  = Start_CYear
    TimTick.CMonth = Start_CMonth
    TimTick.CDay   = Start_CDay
    
    !!! get initial Date
    TimTick.Date     = Math_YMD_To_Date(TimTick.CYear,TimTick.CMonth,TimTick.CDay)
    TimTick.Out_Date = TimTick.Date 
    
    !!!  get Run_AbsNT & Run_AbsNDay if Calender_Run
    call Get_TimeTick_Run_AbsNT_AbsNDay (TimTick)
    
    !!! very important to determine writing restart files frequency
    TimTick.Res_FLAG_Count = 0; !!! be careful
    TimTick.Res_FLAG_Total_Count = 0;
    
    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	write(*,*),'------------------------------------------------------------------------------'
    
    IF( TimTick.Run_Cycle_ID == 1 ) THEN
    
        !!!! CHAR(0) is necessary for C string termination
        Tsub_modloc = trim(Tsub_model_loc)//CHAR(0)
        TimTick.Tsub_model = tsub_new(Tsub_modloc, TimTick.use_gpu) ! initialize the tsub model
        print *, "Tsub Torch Model Initialized"
    
        Tau_modloc = trim(Tau_model_loc)//CHAR(0)
        TimTick.Tau_model = tau_new(Tau_modloc, TimTick.use_gpu)    ! intialize the tau model
	    print *, "Tau Torch Model Initialized"
        
        
    END IF !!! TimTick.Run_Cycle_ID == 1
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    write(*,*),'------------------------------------------------------------------------------'
    write(*,*),'TimeTick Initialization'
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

END SUBROUTINE Alloc_TimeTick 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Alloc_AllFrc_OceanStat (GrdInfo, AllFrc_OcnStat)

    TYPE(GridInfo) :: GrdInfo
    TYPE(AllFrc_OceanStat) :: AllFrc_OcnStat
    
    integer*8 :: NX,NY,NZ,NMode
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    NZ = GrdInfo.NZ
    NMode = GrdInfo.NMode
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate( AllFrc_OcnStat.AllFrc_Co_Pres(NX+2,NY+2,NMode));
    AllFrc_OcnStat.AllFrc_Co_Pres = 0.0d0;
    
    allocate( AllFrc_OcnStat.AllFrc_Co_Ucur(NX+1,NY,NMode));
    AllFrc_OcnStat.AllFrc_Co_Ucur = 0.0d0;
    
    allocate( AllFrc_OcnStat.AllFrc_Co_Vcur(NX,NY+1,NMode));
    AllFrc_OcnStat.AllFrc_Co_Vcur = 0.0d0;
    
    !!! SSTA tendency equation feedbacks
    !!! UATB, UBTA, UATA
    allocate( AllFrc_OcnStat.UATBX(NX+2,NY+2) )
    AllFrc_OcnStat.UATBX = 0.0d0;
    
    allocate( AllFrc_OcnStat.UBTAX(NX+2,NY+2) )
    AllFrc_OcnStat.UBTAX = 0.0d0;
    
    allocate( AllFrc_OcnStat.UATAX(NX+2,NY+2) )
    AllFrc_OcnStat.UATAX = 0.0d0;
    
    !!! VATB, VBTA, VATA
    allocate( AllFrc_OcnStat.VATBY(NX+2,NY+2) )
    AllFrc_OcnStat.VATBY = 0.0d0;
    
    allocate( AllFrc_OcnStat.VBTAY(NX+2,NY+2) )
    AllFrc_OcnStat.VBTAY = 0.0d0;
    
    allocate( AllFrc_OcnStat.VATAY(NX+2,NY+2) )
    AllFrc_OcnStat.VATAY = 0.0d0;
    
    !!! WATB, WBTA, WATA
    allocate( AllFrc_OcnStat.WATBZ(NX+2,NY+2) )
    AllFrc_OcnStat.WATBZ = 0.0d0;
    
    allocate( AllFrc_OcnStat.WBTAZ(NX+2,NY+2) )
    AllFrc_OcnStat.WBTAZ = 0.0d0;
    
    allocate( AllFrc_OcnStat.WATAZ(NX+2,NY+2) )
    AllFrc_OcnStat.WATAZ = 0.0d0;
    
    allocate( AllFrc_OcnStat.SSTA_RESQ(NX+2,NY+2) )
    AllFrc_OcnStat.SSTA_RESQ = 0.0d0;
    
    
    !!! Simple HeatBDG FeedBack 
    allocate( AllFrc_OcnStat.ZAFK(NX+2,NY+2) ); 
    AllFrc_OcnStat.ZAFK = 0.0d0;
    
    allocate( AllFrc_OcnStat.THFK(NX+2,NY+2) ); 
    AllFrc_OcnStat.THFK = 0.0d0;
    
    allocate( AllFrc_OcnStat.RESFK(NX+2,NY+2) ); 
    AllFrc_OcnStat.RESFK = 0.0d0;
    
    allocate( AllFrc_OcnStat.NUDGE(NX+2,NY+2) ); 
    AllFrc_OcnStat.NUDGE = 0.0d0;
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! Simple HeatBDG FeedBack 
    allocate( AllFrc_OcnStat.MonMe_ZAFK(NX+2,NY+2) );   AllFrc_OcnStat.MonMe_ZAFK = 0.0d0;
    allocate( AllFrc_OcnStat.MonMe_THFK(NX+2,NY+2) );   AllFrc_OcnStat.MonMe_THFK = 0.0d0;
    allocate( AllFrc_OcnStat.MonMe_RESFK(NX+2,NY+2) );  AllFrc_OcnStat.MonMe_RESFK = 0.0d0;
    allocate( AllFrc_OcnStat.MonMe_NUDGE(NX+2,NY+2) );  AllFrc_OcnStat.MonMe_NUDGE = 0.0d0;


END SUBROUTINE Alloc_AllFrc_OceanStat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE  Alloc_SurfaceFlux(GrdInfo, SrfFlux)

    TYPE(GridInfo) :: GrdInfo
    TYPE(SurfaceFlux) :: SrfFlux
    
    integer*8 :: NX,NY,NZ,NMode
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    NZ = GrdInfo.NZ
    NMode = GrdInfo.NMode
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    allocate( SrfFlux.TAUX(NX+1, NY) ); 
    SrfFlux.TAUX = 0.0d0;
    allocate( SrfFlux.TAUY(NX, NY+1) ); 
    SrfFlux.TAUY = 0.0d0;
    
    !!! offline test for tau model
    allocate( SrfFlux.TAUX_Test(NX+1, NY) );
    SrfFlux.TAUX_Test = 0.0d0;
    
    allocate( SrfFlux.TAUY_Test(NX, NY+1) );
    SrfFlux.TAUY_Test = 0.0d0;
    
    SrfFlux.TAUX_N4 = 0.0d0; 
    SrfFlux.TAUX_N34 = 0.0d0; 
    SrfFlux.TAUX_N3 = 0.0d0;
    
    SrfFlux.TAUX_N4_Test = 0.0d0;
    SrfFlux.TAUX_N34_Test = 0.0d0;
    SrfFlux.TAUX_N3_Test = 0.0d0;


    allocate( SrfFlux.SSH(NX+2,NY+2) );
    SrfFlux.SSH = 0.0d0;
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    allocate( SrfFlux.MonMe_TAUX(NX+1,NY) ); SrfFlux.MonMe_TAUX(NX+1,NY) = 0.0d0;
    allocate( SrfFlux.MonMe_TAUY(NX,NY+1) ); SrfFlux.MonMe_TAUY(NX,NY+1) = 0.0d0;
    
    allocate( SrfFlux.MonMe_TAUX_Test(NX+1,NY) ); SrfFlux.MonMe_TAUX_Test(NX+1,NY) = 0.0d0;
    allocate( SrfFlux.MonMe_TAUY_Test(NX,NY+1) ); SrfFlux.MonMe_TAUY_Test(NX,NY+1) = 0.0d0;
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! CNOP MJO Perturbation
    
    allocate( SrfFlux.MJO_PsTauX(NX+1,NY) ); SrfFlux.MJO_PsTauX = 0.0d0;
    allocate( SrfFlux.MJO_PsTauY(NX,NY+1) ); SrfFlux.MJO_PsTauY = 0.0d0;
    
    allocate( SrfFlux.PC_MJO_PsTauXY(SrfFlux.NEEOF) ); SrfFlux.PC_MJO_PsTauXY = 0.0d0;
    allocate( SrfFlux.MJO_KelCoeff(SrfFlux.CFM_NOcn) ); SrfFlux.MJO_KelCoeff  = 0.0d0;
    

END SUBROUTINE Alloc_SurfaceFlux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE  Alloc_SurfaceFlux_Hist(GrdInfo, SrfFlux_Hist,HistNT)

    TYPE(GridInfo) :: GrdInfo
    TYPE(SurfaceFlux_Hist) :: SrfFlux_Hist
    
    integer*8 :: NX,NY,NZ,NMode
    integer*8 :: HistNT
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    NZ = GrdInfo.NZ
    NMode = GrdInfo.NMode
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! usually one-year external data (considering leap year)
    
    allocate( SrfFlux_Hist.TAUX_Hist(NX+1, NY, HistNT) ); 
    allocate( SrfFlux_Hist.TAUY_Hist(NX, NY+1, HistNT) );
    
    SrfFlux_Hist.TAUX_Hist = 0.0d0;
    SrfFlux_Hist.TAUY_Hist = 0.0d0;
    
    allocate( SrfFlux_Hist.Nino4_Hist(HistNT) );
    allocate( SrfFlux_Hist.Nino34_Hist(HistNT) );
    allocate( SrfFlux_Hist.Nino3_Hist(HistNT) );
    
    SrfFlux_Hist.Nino4_Hist = 0.0d0;
    SrfFlux_Hist.Nino34_Hist = 0.0d0;
    SrfFlux_Hist.Nino3_Hist = 0.0d0;

    allocate( SrfFlux_Hist.SSH_Hist(NX+2,NY+2,HistNT) );
    SrfFlux_Hist.SSH_Hist = 0.0d0;


END SUBROUTINE Alloc_SurfaceFlux_Hist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Find Nino4/Nino3p4/Nino3 Zone West/East/North/South BoundIndex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Find_Nino_BoundIndex(Grid) 
    
    TYPE(GridInfo) :: Grid
    integer*8 :: IY,IX
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO IY = 2,Grid.NY+1
        
        !!! North Bound Index
        IF( (Grid.PLat(IY+1) > 5.0d0) .AND. ( Grid.PLat(IY) <= 5.0d0 ) )THEN
                  
            Grid.Nino4_Bd.NInd = IY;
            Grid.Nino3p4_Bd.NInd = IY;
            Grid.Nino3_Bd.NInd = IY;
           
        END IF
        
        !!! Equatorial North Bound Index
        IF( (Grid.PLat(IY+1) > 1.0d0) .AND. ( Grid.PLat(IY) <= 1.0d0 ) )THEN
                  
            Grid.Eq_Nino4_Bd.NInd = IY;
            Grid.Eq_Nino3p4_Bd.NInd = IY;
            Grid.Eq_Nino3_Bd.NInd = IY;
           
        END IF
           
        !!! South Bound Index
        IF( (Grid.PLat(IY) >= -5.0d0) .AND. ( Grid.PLat(IY-1) < -5.0d0 ) )THEN
                  
            Grid.Nino4_Bd.SInd = IY;
            Grid.Nino3p4_Bd.SInd = IY;
            Grid.Nino3_Bd.SInd = IY;
           
        END IF
        
        !!! Equatorial South Bound Index
        IF( (Grid.PLat(IY) >= -1.0d0) .AND. ( Grid.PLat(IY-1) < -1.0d0 ) )THEN
                  
            Grid.Eq_Nino4_Bd.SInd = IY;
            Grid.Eq_Nino3p4_Bd.SInd = IY;
            Grid.Eq_Nino3_Bd.SInd = IY;
           
        END IF
           
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO IX = 2,Grid.NX+1
            
        !!! Nino4 West East Bound Index
        IF( ( Grid.PLon(IX-1)< 160.0d0 ).AND.( Grid.PLon(IX) >= 160.0d0 ) )THEN
            Grid.Nino4_Bd.WInd = IX
        END IF
        IF( ( Grid.PLon(IX)<= 210.0d0 ).AND.( Grid.PLon(IX+1) > 210.0d0 ) )THEN
            Grid.Nino4_Bd.EInd = IX
        END IF
            
        !!! Nino3p4 West East Bound Index
        IF( ( Grid.PLon(IX-1)< 190.0d0 ).AND.( Grid.PLon(IX) >= 190.0d0 ) )THEN
            Grid.Nino3p4_Bd.WInd = IX 
        END IF
        IF( ( Grid.PLon(IX)<= 240.0d0 ).AND.( Grid.PLon(IX+1) > 240.0d0 ) )THEN
            Grid.Nino3p4_Bd.EInd = IX
        END IF
            
        !!! Nino3 West East Bound Index
        IF( ( Grid.PLon(IX-1)< 210.0d0 ).AND.( Grid.PLon(IX) >= 210.0d0 ) )THEN
            Grid.Nino3_Bd.WInd = IX 
        END IF
        IF( ( Grid.PLon(IX)<= 270.0d0 ).AND.( Grid.PLon(IX+1) > 270.0d0 ) )THEN
            Grid.Nino3_Bd.EInd = IX
        END IF
            
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!! Equatorial Boundary West East Index
    
    Grid.Eq_Nino4_Bd.WInd = Grid.Nino4_Bd.WInd;
    Grid.Eq_Nino4_Bd.EInd = Grid.Nino4_Bd.EInd;
    
    Grid.Eq_Nino3p4_Bd.WInd = Grid.Nino3p4_Bd.WInd;
    Grid.Eq_Nino3p4_Bd.EInd = Grid.Nino3p4_Bd.EInd;
    
    Grid.Eq_Nino3_Bd.WInd = Grid.Nino3_Bd.WInd;
    Grid.Eq_Nino3_Bd.EInd = Grid.Nino3_Bd.EInd;
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! Nino_PointNum
    
    Grid.Nino4_PointNum = &
    ( Grid.Nino4_Bd.NInd - Grid.Nino4_Bd.SInd +1 ) * ( Grid.Nino4_Bd.EInd - Grid.Nino4_Bd.WInd +1 );
        
    Grid.Nino3p4_PointNum = &
    ( Grid.Nino3p4_Bd.NInd - Grid.Nino3p4_Bd.SInd +1 ) * ( Grid.Nino3p4_Bd.EInd - Grid.Nino3p4_Bd.WInd +1 );
        
    Grid.Nino3_PointNum = &
    ( Grid.Nino3_Bd.NInd - Grid.Nino3_Bd.SInd +1 ) * ( Grid.Nino3_Bd.EInd - Grid.Nino3_Bd.WInd +1 );
    
    !!! Eq_Nino_PointNum
    
    Grid.Eq_Nino4_PointNum = &
    ( Grid.Eq_Nino4_Bd.NInd - Grid.Eq_Nino4_Bd.SInd +1 ) * ( Grid.Eq_Nino4_Bd.EInd - Grid.Eq_Nino4_Bd.WInd +1 );
        
    Grid.Eq_Nino3p4_PointNum = &
    ( Grid.Eq_Nino3p4_Bd.NInd - Grid.Eq_Nino3p4_Bd.SInd +1 ) * ( Grid.Eq_Nino3p4_Bd.EInd - Grid.Eq_Nino3p4_Bd.WInd +1 );
        
    Grid.Eq_Nino3_PointNum = &
    ( Grid.Eq_Nino3_Bd.NInd - Grid.Eq_Nino3_Bd.SInd +1 ) * ( Grid.Eq_Nino3_Bd.EInd - Grid.Eq_Nino3_Bd.WInd +1 );
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
END SUBROUTINE Find_Nino_BoundIndex 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Get_TimeTick_Run_AbsNT_AbsNDay (TimTick)
    !!! get TimTick.AbsNT
    TYPE(TimeTick) :: TimTick
    integer*8 :: DYear, IYear, DMonth
    integer*8 :: NDay
    
    integer*8 :: Start_CYear_Day_MonSum (12), End_CYear_Day_MonSum(12)
    
    DYear = End_CYear - Start_CYear
    DMonth = End_CMonth - Start_CMonth
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (Math_IsLeapYear (Start_CYear) == 0)THEN
    
        Start_CYear_Day_MonSum = Norm_DayMonSum
    
    ELSEIF (Math_IsLeapYear (Start_CYear) == 1)THEN
    
        Start_CYear_Day_MonSum = Leap_DayMonSum
        
    ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (Math_IsLeapYear (End_CYear) == 0)THEN
    
        End_CYear_Day_MonSum = Norm_DayMonSum
    
    ELSEIF (Math_IsLeapYear (Start_CYear) == 1)THEN
    
        End_CYear_Day_MonSum = Leap_DayMonSum
        
    ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF( Leap_Calendar_OnOff == 0 ) THEN
    
        Start_CYear_Day_MonSum = Norm_DayMonSum
        End_CYear_Day_MonSum = Norm_DayMonSum
    
    ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    NDay = 0;
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (Dyear >= 1) THEN !!! run acorss one year
        !!! calculate residual days for Start_CYear
        NDay = NDay + (365 + Leap_Calendar_OnOff * Math_IsLeapYear(Start_CYear) ) &
                    - (Start_CYear_Day_MonSum (Start_CMonth) + (Start_CDay-1) );
    
        !!! calculate days for End_CYear
        NDay = NDay + (End_CYear_Day_MonSum (End_CMonth) + (End_CDay) );
        
    ELSEIF ( DYear == 0 ) THEN !!! run in a same year
    
        NDay = (End_CYear_Day_MonSum (End_CMonth) + (End_CDay) ) - &
               (Start_CYear_Day_MonSum (Start_CMonth) + (Start_CDay-1) );
    
    ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! run acorss more than two year
    IF( DYear >= 2 )THEN 
       DO IYear = (Start_CYear+1),(End_CYear-1)
          NDay = NDay + 365 + Leap_Calendar_OnOff * Math_IsLeapYear(IYear);
       END DO
    ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TimTick.Run_AbsNDay = NDay;
    TimTick.Run_AbsNT = NDay * (Day_NT);

END SUBROUTINE Get_TimeTick_Run_AbsNT_AbsNDay
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Update_TimeTick (TimTick)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! can you calculate CYear,CMonth,CDay 
!!!          based on Run_AbsIT, Start_CYear, Start_CMonth, Start_CDay ?
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TYPE(TimeTick) :: TimTick
    
    integer*8 :: IDay,New_IDay,Old_IDay
    integer*8 :: DayMon(12)
    
    !!! reached the end of day/month/year
    !integer*8 :: Reach_CDay_FLAG = 0, Reach_CMonth_FLAG = 0, Reach_CYear_FLAG = 0
    integer*8 :: Reach_FLAG = 0
    
    !!!! arrived at the start of day/month/year 
    !integer*8 :: Arrive_CDay_FLAG = 0, Arrive_CMonth_FLAG=0, Arrive_CYear_FLAG = 0
    !integer*8 :: Arrive_FLAG = 0
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    IDay = ceiling(dble(TimTick.Run_AbsIT)/dble(Day_NT));
    TimTick.Run_AbsIDay = IDay;
    
    New_IDay = ceiling((dble(TimTick.Run_AbsIT+1))/dble(Day_NT));
    Old_IDay = ceiling((dble(TimTick.Run_AbsIT-1))/dble(Day_NT));
    IF (TimTick.Run_AbsIT == 1 )THEN
        Old_IDay = 1
    END IF !!! for safety, so (IDay - Old_IDay).EQ.0
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (Math_IsLeapYear (TimTick.CYear) == 0)THEN
       DayMon = Norm_DayMon 
    ELSEIF (Math_IsLeapYear (TimTick.CYear) == 1)THEN
       DayMon = Leap_DayMon
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF( Leap_Calendar_OnOff == 0 ) THEN
        DayMon = Norm_DayMon
    ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( (IDay - Old_IDay).NE.0 )THEN 
    !!! day increment
         TimTick.CDay = TimTick.CDay + 1;
         !!! whether surpass a month
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF ( TimTick.CDay > DayMon(TimTick.CMonth) ) THEN
              TimTick.CDay = 1;
              TimTick.CMonth = TimTick.CMonth + 1;
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              IF (TimTick.CMonth > 12 )THEN
                 TimTick.CMonth = 1;
                 TimTick.CYear = TimTick.CYear + 1;
              END IF !!! cross year
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         END IF !!! cross month
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END IF !!! cross day
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! update Date
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TimTick.Date = Math_YMD_To_Date(TimTick.CYear,TimTick.CMonth,TimTick.CDay)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! update Run_AbsIMon 
    IF( ( TimTick.CYear - Start_CYear ) == 0 ) THEN 
        TimTick.Run_AbsIMon = (TimTick.CMonth - Start_CMonth + 1);
    ELSEIF( ( TimTick.CYear - Start_CYear ) == 1 ) THEN 
        TimTick.Run_AbsIMon = (12 - Start_CMonth + 1) + TimTick.CMonth;
    ELSE
        TimTick.Run_AbsIMon = &
        (12 - Start_CMonth + 1) + 12*(TimTick.CYear - Start_CYear - 1) + TimTick.CMonth;
    ENDIF
    
    !!! set Reach_CDay_FLAG
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( mod(TimTick.Run_AbsIT,Day_NT)==0 )THEN
       TimTick.Reach_CDay_FLAG = 1;
    ELSE
       TimTick.Reach_CDay_FLAG = 0;
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! set Reach_CMonth_FLAG
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( (TimTick.CDay == DayMon(TimTick.CMonth)) .AND. (mod(TimTick.Run_AbsIT,Day_NT)==0) )THEN
       TimTick.Reach_CMonth_FLAG = 1;
    ELSE
       TimTick.Reach_CMonth_FLAG = 0;
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! set Reach_CYear_FLAG
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( (TimTick.CMonth == 12).AND. (TimTick.CDay == 31) .AND. (mod(TimTick.Run_AbsIT,Day_NT)==0))THEN
       TimTick.Reach_CYear_FLAG = 1;
    ELSE
       TimTick.Reach_CYear_FLAG = 0;
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!! set Arrive_CDay_FLAG
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( mod(TimTick.Run_AbsIT,Day_NT)==1 )THEN
       TimTick.Arrive_CDay_FLAG = 1;
    ELSE
       TimTick.Arrive_CDay_FLAG = 0;
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! set Arrive_CMonth_FLAG
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( (TimTick.CDay == 1).AND. (mod(TimTick.Run_AbsIT,Day_NT)==1) )THEN
       TimTick.Arrive_CMonth_FLAG = 1;
    ELSE
       TimTick.Arrive_CMonth_FLAG = 0;
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! set Arrive_CYear_FLAG
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( (TimTick.CMonth == 1).AND. (TimTick.CDay == 1) .AND. (mod(TimTick.Run_AbsIT,Day_NT)==1))THEN
       TimTick.Arrive_CYear_FLAG = 1;
    ELSE
       TimTick.Arrive_CYear_FLAG = 0;
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! we need these FLAGs to generate Restart File
    IF( (ReStart_DT_Unit == 'D') )THEN
        
        IF (TimTick.Reach_CDay_FLAG == 1 ) THEN
        
            TimTick.Res_FLAG_Count = TimTick.Res_FLAG_Count + 1
            TimTick.Res_FLAG_Total_Count = TimTick.Res_FLAG_Total_Count + 1
      
        ENDIF
        
    ELSEIF( (ReStart_DT_Unit == 'M') )THEN
    
        IF (TimTick.Reach_CMonth_FLAG == 1 ) THEN
            
            TimTick.Res_FLAG_Count = TimTick.Res_FLAG_Count + 1
            TimTick.Res_FLAG_Total_Count = TimTick.Res_FLAG_Total_Count + 1
        
        ENDIF
        
    ELSEIF( (ReStart_DT_Unit == 'Y') )THEN
    
        IF (TimTick.Reach_CYear_FLAG == 1 ) THEN
            
            TimTick.Res_FLAG_Count = TimTick.Res_FLAG_Count + 1
            TimTick.Res_FLAG_Total_Count = TimTick.Res_FLAG_Total_Count + 1
      
        ENDIF
        
    END IF
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! get YrDaySum
    IF ( Leap_Calendar_OnOff == 1) THEN
        IF ( Math_IsLeapYear(TimTick.CYear) == 1) THEN
        TimTick.YrDaySum = (Leap_DayMonSum (TimTick.CMonth) + (TimTick.CDay) );
        ELSE
        TimTick.YrDaySum = (Norm_DayMonSum (TimTick.CMonth) + (TimTick.CDay) );
        END IF
    ELSE
        TimTick.YrDaySum = (Norm_DayMonSum (TimTick.CMonth) + (TimTick.CDay) );
    END IF
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Show TimeTick on screen
    IF(Show_TimeTick_DT_Unit == 'D')THEN
    
       Reach_FLAG = TimTick.Reach_CDay_FLAG
       
    ELSEIF(Show_TimeTick_DT_Unit == 'M')THEN
    
       Reach_FLAG = TimTick.Reach_CMonth_FLAG
       
    ELSEIF(Show_TimeTick_DT_Unit == 'Y')THEN
    
       Reach_FLAG = TimTick.Reach_CYear_FLAG
       
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (Reach_FLAG == 1) THEN
    
    write(6,200) TimTick.Cyear,TimTick.Cmonth,TimTick.CDay, TimTick.Run_AbsIDay, TimTick.Run_AbsIMon
    200 format(" Cyear = ",I4,"  Cmonth = ",I2,"  CDay = ",I2, " Run_AbsIDay =",I8, "  Run_AbsIMon =",I6)
    
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
END SUBROUTINE Update_TimeTick
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Read_Grid_ClimStat (GrdPath, ClimPath,FrcPath, GrdInfo, Clim_OcnStat)

    character(len=*) :: GrdPath, ClimPath, FrcPath
    TYPE(GridInfo) :: GrdInfo
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    
    real*8 :: dNX,dNY,dNZ,dNMode !!!! double to integer
    integer*8 :: NX,NY,NZ,NMode
    integer*8 :: IX,IY,IZ,IM,IT
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! tmp Wsrf,TAUX,TAUY to compute model climtalogical Wsrf
    real*8,allocatable :: tmp_Wsrf(:,:),tmp_TAUX(:,:),tmp_TAUY(:,:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*),'---------------------------------------------------------------------------'
    write(*,*),'Reading GridInfo ...'
    write(*,*),''
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call NC_Read_One_IT (GrdPath,'PGrid.nc','NX',dNX,1,'A');
    GrdInfo.NX = int8(dNX);
    NX = int8(dNX);
    
    call NC_Read_One_IT (GrdPath,'PGrid.nc','NY',dNY,1,'A');
    GrdInfo.NY = int8(dNY);
    NY = int8(dNY);
    
    call NC_Read_One_IT (GrdPath,'PGrid.nc','NZ',dNZ,1,'A');
    GrdInfo.NZ = int8(dNZ);
    NZ = int8(dNZ);
    
    call NC_Read_One_IT (GrdPath,'PGrid.nc','NMode',dNMode,1,'A');
    GrdInfo.NMode = int8(dNMode);
    NMode = int8(dNMode);
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    allocate (GrdInfo.IMode(NMode));
    DO IM = 1,NMode
       GrdInfo.IMode(IM) = dble(IM)
    END DO
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call NC_Read_Mul_IT (GrdPath,'PGrid.nc','Depth',GrdInfo.Depth,0,0,'A');
    
    call NC_Read_Mul_IT (GrdPath,'PGrid.nc','PLon',GrdInfo.PLon,0,0,'A');
    call NC_Read_Mul_IT (GrdPath,'PGrid.nc','PLat',GrdInfo.PLat,0,0,'A');
    call NC_Read_Mul_IT (GrdPath,'UGrid.nc','ULon',GrdInfo.ULon,0,0,'A');
    call NC_Read_Mul_IT (GrdPath,'UGrid.nc','ULat',GrdInfo.ULat,0,0,'A');
    call NC_Read_Mul_IT (GrdPath,'VGrid.nc','VLon',GrdInfo.VLon,0,0,'A');
    call NC_Read_Mul_IT (GrdPath,'VGrid.nc','VLat',GrdInfo.VLat,0,0,'A');
    
    call NC_Read_Mul_IT (GrdPath,'PGrid.nc','PReX',GrdInfo.PReX,0,0,'A');
    call NC_Read_Mul_IT (GrdPath,'PGrid.nc','PReY',GrdInfo.PReY,0,0,'A');
    call NC_Read_Mul_IT (GrdPath,'UGrid.nc','UReX',GrdInfo.UReX,0,0,'A');
    call NC_Read_Mul_IT (GrdPath,'UGrid.nc','UReY',GrdInfo.UReY,0,0,'A');
    call NC_Read_Mul_IT (GrdPath,'VGrid.nc','VReX',GrdInfo.VReX,0,0,'A');
    call NC_Read_Mul_IT (GrdPath,'VGrid.nc','VReY',GrdInfo.VReY,0,0,'A');
    
    call NC_Read_Mul_IT (GrdPath,'PGrid.nc','PMask',GrdInfo.PMask,0,0,'A');
    call NC_Read_Mul_IT (GrdPath,'UGrid.nc','UMask',GrdInfo.UMask,0,0,'A');
    call NC_Read_Mul_IT (GrdPath,'VGrid.nc','VMask',GrdInfo.VMask,0,0,'A');
    
    !!!! e-folding Pmask to preprocess input for TAUXY_AI
    !!!! call NC_Read_Mul_IT (GrdPath,'PGrid.nc','Tau_PMask',GrdInfo.Tau_PMask,0,0,'A');
    allocate(GrdInfo.Tau_PMask(size(GrdInfo.PMask,1),size(GrdInfo.PMask,2)))
    GrdInfo.Tau_PMask = 0.0d0;
    
    DO IX = 1,size(GrdInfo.PMask,1)
    
       GrdInfo.Tau_PMask(IX,1:NY+2) = GrdInfo.PMask(IX,1:NY+2) * &
       exp( (-EarthBeta/(2.0d0 * EFold_C0)) * (GrdInfo.PReY)**(2.0d0) );
    
    END DO
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! get ULNZ & ULDepth 
    DO IZ = 1,(GrdInfo.NZ-1)
    IF ( (GrdInfo.Depth(IZ) <= Ref_UppLayer_Depth).AND.( GrdInfo.Depth(IZ+1) > Ref_UppLayer_Depth) )THEN
        GrdInfo.ULNZ = IZ;
    END IF
    END DO
    
    allocate(GrdInfo.ULDepth(GrdInfo.ULNZ));
    GrdInfo.ULDepth = GrdInfo.Depth(1:GrdInfo.ULNZ);
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    GrdInfo.ReDX = abs( GrdInfo.PReX(2)-GrdInfo.PReX(1) );
    GrdInfo.ReDY = abs( GrdInfo.PReY(2)-GrdInfo.PReY(1) );
    
    GrdInfo.DLon = abs( GrdInfo.PLon(2)-GrdInfo.PLon(1) );
    GrdInfo.DLat = abs( GrdInfo.PLat(2)-GrdInfo.PLat(1) );
    GrdInfo.ReDT = Ocn_NumDT; !!! be careful
    
    GrdInfo.ReDZ = abs( GrdInfo.Depth(2)-GrdInfo.Depth(1) );
    
    !!!! Get Nino BoundIndex
    call Find_Nino_BoundIndex (GrdInfo);
    
    !!!! Set Sponge Layers
    allocate(GrdInfo.Sponge_Alp(GrdInfo.Sponge_Num));
    DO IY = 1,GrdInfo.Sponge_Num
       GrdInfo.Sponge_Alp(IY) = ( 1.0d0 - dble(IY)/dble(GrdInfo.Sponge_Num) ) ** dble(GrdInfo.Sponge_Pow)
    END DO
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*),'GridInfo Reading Completed'
    write(*,*),'-----------------------------------------------------------------------------'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*),'-----------------------------------------------------------------------------'
    write(*,*),'Reading Clim Eigen Baroclinic Mode ...'
    write(*,*),''
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call NC_Read_Mul_IT (ClimPath,EigMode_FName,'EigMode',Clim_OcnStat.EigMode,0,0,'A');
    call NC_Read_Mul_IT (ClimPath,EigMode_FName,'EigC',Clim_OcnStat.EigC,0,0,'A');
    call NC_Read_Mul_IT (ClimPath,EigMode_FName,'Co_Wstr',Clim_OcnStat.Co_Wstr,0,0,'A');
    
    call NC_Read_Mul_IT (ClimPath,EigMode_FName,'IZ_EigMode',Clim_OcnStat.IZ_EigMode,0,0,'A');
    call NC_Read_Mul_IT (ClimPath,EigMode_FName,'DZ_EigMode',Clim_OcnStat.DZ_EigMode,0,0,'A');
    
    !!! get ClimNT
    GrdInfo.ClimNT = 12; 
    
    !!! get (EigC).^2
    allocate(Clim_OcnStat.EigC_Pow2(NX+2,NMode));
    Clim_OcnStat.EigC_Pow2 = Clim_OcnStat.EigC * Clim_OcnStat.EigC;
    
    !!! write(*,*),'Co_Wstr =',Clim_OcnStat.Co_Wstr
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*),'Clim Baroclinic Mode Reading Completed...'
    write(*,*),'-----------------------------------------------------------------------------'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*),'-----------------------------------------------------------------------------'
    write(*,*),'Reading Clim PsTauXY OceanStat and Computing Clim TauXY ...'
    write(*,*),''
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    call NC_Read_Mul_IT (ClimPath,Clim_PsTauX_FName,'PsTauX',Clim_OcnStat.PsTauXBar,0,0,'A');
    call NC_Read_Mul_IT (ClimPath,Clim_PsTauY_FName,'PsTauY',Clim_OcnStat.PsTauYBar,0,0,'A');
    
    call NC_Read_Mul_IT (ClimPath,CMIP_Clim_TauX_FName,'CMIP_TauX',Clim_OcnStat.CMIP_TauXBar,0,0,'A');
    call NC_Read_Mul_IT (ClimPath,CMIP_Clim_TauY_FName,'CMIP_TauY',Clim_OcnStat.CMIP_TauYBar,0,0,'A');
    
    !!! compute clim TauX and TauY
    !allocate( Clim_OcnStat.TauXBar (NX+1,NY,GrdInfo.ClimNT) );
    !Clim_OcnStat.TauXBar = CoeffDrag * RhoAtm * Clim_OcnStat.PsTauXBar;
    
    !allocate( Clim_OcnStat.TauYBar (NX,NY+1,GrdInfo.ClimNT) );
    !Clim_OcnStat.TauYBar = CoeffDrag * RhoAtm * Clim_OcnStat.PsTauYBar;

    call NC_Read_Mul_IT (ClimPath,Clim_TauX_FName,'UFLX',Clim_OcnStat.TauXBar,0,0,'A');
    call NC_Read_Mul_IT (ClimPath,Clim_TauY_FName,'VFLX',Clim_OcnStat.TauYBar,0,0,'A');
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*),'Clim PsTauXY Reading and TauXY Computing Completed'
    write(*,*),'-----------------------------------------------------------------------------'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*),'-----------------------------------------------------------------------------'
    write(*,*),'Reading Clim Usrf, Vsrf, Wsrf, SSH, SST, TSUB, MLD, TCD ...'
    write(*,*),''
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    call NC_Read_Mul_IT (ClimPath,Clim_OceanStat_FName,'UsrfBar',Clim_OcnStat.UsrfBar,0,0,'A');
    call NC_Read_Mul_IT (ClimPath,Clim_OceanStat_FName,'VsrfBar',Clim_OcnStat.VsrfBar,0,0,'A');
    
    call NC_Read_Mul_IT (ClimPath,Clim_OceanStat_FName,'WsrfBar',Clim_OcnStat.WsrfBar,0,0,'A');
    call NC_Read_Mul_IT (ClimPath,Clim_OceanStat_FName,'MLDBar',Clim_OcnStat.MLDBar,0,0,'A');
    call NC_Read_Mul_IT (ClimPath,Clim_OceanStat_FName,'TCDBar',Clim_OcnStat.TCDBar,0,0,'A');
    
    call NC_Read_Mul_IT (ClimPath,Clim_OceanStat_FName,'SSHBar',Clim_OcnStat.SSHBar,0,0,'A');
    call NC_Read_Mul_IT (ClimPath,Clim_OceanStat_FName,'SSTBar',Clim_OcnStat.SSTBar,0,0,'A');
    call NC_Read_Mul_IT (ClimPath,Clim_OceanStat_FName,'TSUBBar',Clim_OcnStat.TSUBBar,0,0,'A');
    
    call NC_Read_Mul_IT (ClimPath,Clim_OceanStat_FName,'CMIP_SSHBar',Clim_OcnStat.CMIP_SSHBar,0,0,'A');
    call NC_Read_Mul_IT (ClimPath,Clim_OceanStat_FName,'CMIP_SSTBar',Clim_OcnStat.CMIP_SSTBar,0,0,'A');
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! MASK for safety
    !DO IT = 1, GrdInfo.ClimNT
    !
    !Clim_OcnStat.UsrfBar(1:NX+1,1:NY,IT) = Clim_OcnStat.UsrfBar(1:NX+1,1:NY,IT) * GrdInfo.UMask
    !Clim_OcnStat.VsrfBar(1:NX,1:NY+1,IT) = Clim_OcnStat.VsrfBar(1:NX,1:NY+1,IT) * GrdInfo.VMask
    !Clim_OcnStat.WsrfBar(1:NX+2,1:NY+2,IT)  = Clim_OcnStat.WsrfBar(1:NX+2,1:NY+2,IT) * GrdInfo.PMask
    !
    !Clim_OcnStat.SSHBar (1:NX+2,1:NY+2,IT)  = Clim_OcnStat.SSHBar (1:NX+2,1:NY+2,IT) * GrdInfo.PMask
    !Clim_OcnStat.SSTBar (1:NX+2,1:NY+2,IT)  = Clim_OcnStat.SSTBar (1:NX+2,1:NY+2,IT) * GrdInfo.PMask
    !
    !!!! MLD & TCD Bar may work for WATBZ and WBTAZ
    !!!! Clim_OcnStat.MLDBar (1:NX+2,1:NY+2,IT)  = Clim_OcnStat.MLDBar (1:NX+2,1:NY+2,IT) * GrdInfo.PMask
    !!!! Clim_OcnStat.TCDBar (1:NX+2,1:NY+2,IT)  = Clim_OcnStat.TCDBar (1:NX+2,1:NY+2,IT) * GrdInfo.PMask
    !
    !END DO 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*),'Clim Usrf, Vsrf, Wsrf, SSH, SST, TSUB, MLD, TCD Reading completed'
    write(*,*),'-----------------------------------------------------------------------------'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*),'-----------------------------------------------------------------------------'
    write(*,*),'Computing Clim Wsrf according to TauXBar, TauYBar'
    write(*,*),' '
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate( tmp_Wsrf(NX+2,NY+2), tmp_TAUX(NX+1,NY), tmp_TAUY(NX,NY+1) );
    tmp_Wsrf = 0.0d0;
    tmp_TAUX = 0.0d0;
    tmp_TAUY = 0.0d0;
    
    DO IT = 1, GrdInfo.ClimNT
    
        !tmp_TAUX = CoeffDrag * RhoAtm * Clim_OcnStat.PsTauXBar(1:NX+1,1:NY,IT)
        !tmp_TAUY = CoeffDrag * RhoAtm * Clim_OcnStat.PsTauYBar(1:NX,1:NY+1,IT)

        tmp_TAUX = Clim_OcnStat.TAUXBar(1:NX+1,1:NY,IT)
        tmp_TAUY = Clim_OcnStat.TAUYBar(1:NX,1:NY+1,IT)
        
        !!! Ekman Layer dynamics 
        call Get_Ekman_Wsrf(GrdInfo, tmp_Wsrf, tmp_TAUX, tmp_TAUY)
        
        Clim_OcnStat.WsrfBar(1:NX+2,1:NY+2,IT) = tmp_Wsrf;
    
    END DO
    
    deallocate( tmp_Wsrf, tmp_TAUX, tmp_TAUY );
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*),'Clim Wsrf according to TauXBar, TauYBar computing completed'
    write(*,*),'-----------------------------------------------------------------------------'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !write(*,*),'-----------------------------------------------------------------------------'
    !write(*,*),'Reading EOF Linear Tau Model Parameters ...'
    !write(*,*),''
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !call NC_Read_Mul_IT &
    !     (ClimPath,EOF_Linear_Tau_Model_FName,'SSH_SST_scale',Clim_OcnStat.SSH_SST_scale,0,0,'A');
    !
    !call NC_Read_Mul_IT &
    !     (ClimPath,EOF_Linear_Tau_Model_FName,'SSH_SST_EOF',Clim_OcnStat.SSH_SST_EOF,0,0,'A');
    !
    !call NC_Read_Mul_IT &
    !     (ClimPath,EOF_Linear_Tau_Model_FName,'TAUX_TAUY_EOF',Clim_OcnStat.TAUX_TAUY_EOF,0,0,'A');
    !
    !call NC_Read_Mul_IT &
    !     (ClimPath,EOF_Linear_Tau_Model_FName,'Tau_AMat_12Mon',Clim_OcnStat.Tau_AMat_12Mon,0,0,'A');
    !
    !call NC_Read_Mul_IT &
    !     (ClimPath,EOF_Linear_Tau_Model_FName,'Tau_BVec_12Mon',Clim_OcnStat.Tau_BVec_12Mon,0,0,'A');
    !
    !Clim_OcnStat.Tau_EOF_NMode = size(Clim_OcnStat.Tau_AMat_12Mon,1);
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !write(*,*),'Tau_EOF_NMode        = '//trim(Num2Str(Clim_OcnStat.Tau_EOF_NMode))
    !write(*,*),'Tau_EOF_Trunc_NMode  = '//trim(Num2Str(Tau_EOF_Trunc_NMode))
    !write(*,*),'Tau_PC1_Coeff        = '//trim(Num2Str(Clim_OcnStat.Tau_AMat_12Mon(1,1,1)))
    !write(*,*),''
    !write(*,*),'EOF Linear Tau Model Parameters Reading completed'
    !write(*,*),'-----------------------------------------------------------------------------'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*),'-----------------------------------------------------------------------------'
    write(*,*),'Reading Regression Linear Tau Model Parameters ...'
    write(*,*),''
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!! TAUX
    
    call NC_Read_Mul_IT &
         (ClimPath,Reg_Linear_Tau_Model_FName,'TAUX_Reg_N4',Clim_OcnStat.TAUX_Reg_N4,0,0,'A');
    
    call NC_Read_Mul_IT &
         (ClimPath,Reg_Linear_Tau_Model_FName,'TAUX_Reg_N3',Clim_OcnStat.TAUX_Reg_N3,0,0,'A');
    
    call NC_Read_Mul_IT &
         (ClimPath,Reg_Linear_Tau_Model_FName,'TAUX_Reg_Const',Clim_OcnStat.TAUX_Reg_Const,0,0,'A');
    
    !!! TAUY
    call NC_Read_Mul_IT &
         (ClimPath,Reg_Linear_Tau_Model_FName,'TAUY_Reg_N4',Clim_OcnStat.TAUY_Reg_N4,0,0,'A');
    
    call NC_Read_Mul_IT &
         (ClimPath,Reg_Linear_Tau_Model_FName,'TAUY_Reg_N3',Clim_OcnStat.TAUY_Reg_N3,0,0,'A');
    
    call NC_Read_Mul_IT &
         (ClimPath,Reg_Linear_Tau_Model_FName,'TAUY_Reg_Const',Clim_OcnStat.TAUY_Reg_Const,0,0,'A');
    
    !!! SSTA
    call NC_Read_Mul_IT &
         (ClimPath,Reg_Linear_Tau_Model_FName,'SSTA_Reg_N4',Clim_OcnStat.SSTA_Reg_N4,0,0,'A');
    
    call NC_Read_Mul_IT &
         (ClimPath,Reg_Linear_Tau_Model_FName,'SSTA_Reg_N3',Clim_OcnStat.SSTA_Reg_N3,0,0,'A');
    
    call NC_Read_Mul_IT &
         (ClimPath,Reg_Linear_Tau_Model_FName,'SSTA_Reg_Const',Clim_OcnStat.SSTA_Reg_Const,0,0,'A');
    
    !!! Skewed SSTA
    call NC_Read_Mul_IT &
         (ClimPath,Reg_Linear_Tau_Model_FName,'Skew_SSTA_Reg_N4',Clim_OcnStat.Skew_SSTA_Reg_N4,0,0,'A');
    
    call NC_Read_Mul_IT &
         (ClimPath,Reg_Linear_Tau_Model_FName,'Skew_SSTA_Reg_N3',Clim_OcnStat.Skew_SSTA_Reg_N3,0,0,'A');
    
    call NC_Read_Mul_IT &
         (ClimPath,Reg_Linear_Tau_Model_FName,'Skew_SSTA_Reg_Const',Clim_OcnStat.Skew_SSTA_Reg_Const,0,0,'A');
    
    
    !!! Skewed param
    call NC_Read_One_IT &
         (ClimPath,Reg_Linear_Tau_Model_FName,'N4_Skew_A',Clim_OcnStat.N4_Skew_A,1,'A');
    call NC_Read_One_IT &
         (ClimPath,Reg_Linear_Tau_Model_FName,'N4_Skew_B',Clim_OcnStat.N4_Skew_B,1,'A');
    call NC_Read_One_IT &
         (ClimPath,Reg_Linear_Tau_Model_FName,'N4_Skew_C',Clim_OcnStat.N4_Skew_C,1,'A');
    call NC_Read_One_IT &
         (ClimPath,Reg_Linear_Tau_Model_FName,'N4_Skew_M',Clim_OcnStat.N4_Skew_M,1,'A');
    
    call NC_Read_One_IT &
         (ClimPath,Reg_Linear_Tau_Model_FName,'N3_Skew_A',Clim_OcnStat.N3_Skew_A,1,'A');
    call NC_Read_One_IT &
         (ClimPath,Reg_Linear_Tau_Model_FName,'N3_Skew_B',Clim_OcnStat.N3_Skew_B,1,'A');
    call NC_Read_One_IT &
         (ClimPath,Reg_Linear_Tau_Model_FName,'N3_Skew_C',Clim_OcnStat.N3_Skew_C,1,'A');
    call NC_Read_One_IT &
         (ClimPath,Reg_Linear_Tau_Model_FName,'N3_Skew_M',Clim_OcnStat.N3_Skew_M,1,'A');

    
    !!!! Std_SSTA
    !call NC_Read_Mul_IT &
    !     (ClimPath,Reg_Linear_Tau_Model_FName,'Std_SSTA',Clim_OcnStat.Std_SSTA,0,0,'A');

    !!! THFK_Mask
    call NC_Read_Mul_IT &
        (ClimPath,Reg_Linear_Tau_Model_FName,'THFK_Mask',Clim_OcnStat.THFK_Mask,0,0,'A'); 
    
    !!! TSUB_OutPut_Mask
    call NC_Read_Mul_IT &
        (ClimPath,Reg_Linear_Tau_Model_FName,'TSUB_OutPut_Mask',Clim_OcnStat.TSUB_OutPut_Mask,0,0,'A'); 
    
    !!! EOF TSUB to remove noise of TSUB-AI
    call NC_Read_Mul_IT &
        (ClimPath,Reg_Linear_Tau_Model_FName,'EOF_TSUB',Clim_OcnStat.EOF_TSUB,0,0,'A'); 
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    write(*,*),'size(TAUX_Reg_N4,1) = ',trim(Num2Str(size(Clim_OcnStat.TAUX_Reg_N4,1)));
    write(*,*),'size(TAUX_Reg_N4,2) = ',trim(Num2Str(size(Clim_OcnStat.TAUX_Reg_N4,2)));
    write(*,*),'size(TAUX_Reg_N4,3) = ',trim(Num2Str(size(Clim_OcnStat.TAUX_Reg_N4,3)));
    
    write(*,*),'' 
    write(*,*),'Regression Linear Tau Model Parameters Reading completed'
    write(*,*),'-----------------------------------------------------------------------------'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*),'-----------------------------------------------------------------------------'
    write(*,*),'Reading MJO Forcing EEOF ...'
    write(*,*),''
    
     
    call NC_Read_Mul_IT &
         (FrcPath,MJO_Forcing_EEOF_FName,'EEOF_PsTauX',Clim_OcnStat.EEOF_PsTauX,0,0,'A');
    
    call NC_Read_Mul_IT &
         (FrcPath,MJO_Forcing_EEOF_FName,'EEOF_PsTauY',Clim_OcnStat.EEOF_PsTauY,0,0,'A');
    
    call NC_Read_Mul_IT &
         (FrcPath,MJO_Forcing_EEOF_FName,'PC_PsTauXY',Clim_OcnStat.PC_PsTauXY,0,0,'A');
    
    call NC_Read_Mul_IT &
         (FrcPath,MJO_Forcing_EEOF_FName,'EEOF_KelCoeff',Clim_OcnStat.EEOF_KelCoeff,0,0,'A');
    
    
    write(*,*),'size(EEOF_PsTauX,1) = ',trim(Num2Str(size(Clim_OcnStat.EEOF_PsTauX,1)));
    write(*,*),'size(EEOF_PsTauX,2) = ',trim(Num2Str(size(Clim_OcnStat.EEOF_PsTauX,2)));
    write(*,*),'size(EEOF_PsTauX,3) = ',trim(Num2Str(size(Clim_OcnStat.EEOF_PsTauX,3)));
    write(*,*),'size(EEOF_PsTauX,4) = ',trim(Num2Str(size(Clim_OcnStat.EEOF_PsTauX,4)));
    
    write(*,*),'' 
    write(*,*),'MJO Forcing EEOF Reading completed'
    write(*,*),'-----------------------------------------------------------------------------'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(6,110) MY_OMP_THREADS
    110  format(' MY_OMP_THREADS  = ',I4)
    
    write(6,111) GrdInfo.NX
    111  format(' GridInfo.NX     = ',I4)
    
    write(6,112) GrdInfo.NY
    112  format(' GridInfo.NY     = ',I4)
    
    write(6,113) GrdInfo.NZ
    113  format(' GridInfo.NZ     = ',I4)
    
    write(6,114) GrdInfo.ULNZ
    114  format(' GridInfo.ULNZ   = ',I4)
    
    write(6,115) GrdInfo.ClimNT
    115  format(' GridInfo.ClimNT = ',I4)
    
    write(6,116) GrdInfo.NMode
    116  format(' GridInfo.NMode  = ',I4)
    
    write(6,117) int8(GrdInfo.DLon)
    117  format(' GridInfo.DLon   = ',I4)
    
    write(6,118) int8(GrdInfo.DLat)
    118  format(' GridInfo.DLat   = ',I4)
    
    write(6,119) int8(GrdInfo.ReDZ)
    119  format(' GridInfo.ReDZ   = ',I4)
    
    write(6,120) int8(GrdInfo.Nino4_PointNum)
    120  format(' GridInfo.N4_PN  = ',I5) 
    
    write(6,121) int8(GrdInfo.Nino3p4_PointNum)
    121  format(' GridInfo.N34_PN = ',I5) 
    
    write(6,122) int8(GrdInfo.Nino3_PointNum)
    122  format(' GridInfo.N4_PN  = ',I5) 
    
    write(6,123) int8(GrdInfo.Sponge_Num)
    123  format(' Sponge_Num      = ',I4) 
    
    write(6,124) int8(GrdInfo.Sponge_Pow)
    124  format(' Sponge_Pow      = ',I4) 
    
    write(*,*) 'GridInfo.ReDT   = '//trim(Num2Str(GrdInfo.ReDT/3600.0d0))//' hour'
    write(*,*),'------------------------------------------------------------------------------'
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE Read_Grid_ClimStat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Write_ReStart(ResPath,GrdInfo,OcnStat,TimTick,SrfFlux);

     TYPE (OceanStat) :: OcnStat
     TYPE (GridInfo) :: GrdInfo
     character(len=*) :: ResPath
     TYPE(TimeTick) :: TimTick
     TYPE(SurfaceFlux) :: SrfFlux
     
     character(len=300) :: P_ResName = 'ReStart_PGrid.nc';
     character(len=300) :: U_ResName = 'ReStart_UGrid.nc';
     character(len=300) :: V_ResName = 'ReStart_VGrid.nc';
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     real*8,allocatable :: res_time_tick(:) !!! tmp timetick size == 1
     character(len=100) :: time_unit   !!! time_unit
     
     allocate(res_time_tick(1)); res_time_tick(1) = dble(TimTick.Run_AbsIDay);
     
     !!! to prevent misjudge of CDO (sub half day)
     res_time_tick = res_time_tick - 0.5d0
     
     time_unit = &
     'days since '//trim(Math_YMD_To_Date(Start_CYear,Start_CMonth,Start_CDay))//' 00:00:0.0';
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     write(*,*),'------------------------------------------------------------------------------'
     write(*,*),'Writing ReStart Files ...'

     call  NC_Create_File( ResPath,P_ResName, &
                           GrdInfo.PLon, GrdInfo.PLat, GrdInfo.IMode, res_time_tick, &
                           'lon', 'lat','imode','time',time_unit, 'X', Leap_Calendar_OnOff)
    
     call  NC_Create_File( ResPath,U_ResName, &
                           GrdInfo.ULon, GrdInfo.ULat, GrdInfo.IMode, res_time_tick, &
                           'lon', 'lat', 'imode','time',time_unit, 'X',Leap_Calendar_OnOff)
    
     call  NC_Create_File( ResPath,V_ResName, &
                           GrdInfo.VLon, GrdInfo.VLat, GrdInfo.IMode, res_time_tick, &
                           'lon', 'lat', 'imode','time',time_unit,'X', Leap_Calendar_OnOff)
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     !!! be careful, in restart file, numerical precision should be double
     call NC_Def_Var(ResPath,P_ResName,'Co_Pres','none','lon','lat','imode','time',Double_FLAG=1)
     call NC_Def_Var(ResPath,U_ResName,'Co_Ucur','none','lon','lat','imode','time',Double_FLAG=1)
     call NC_Def_Var(ResPath,V_ResName,'Co_Vcur','none','lon','lat','imode','time',Double_FLAG=1)
     
     call NC_Def_Var(ResPath,P_ResName,'SST','dC','lon','lat','time',Double_FLAG=1)
     call NC_Def_Var(ResPath,P_ResName,'TSUB','dC','lon','lat','time',Double_FLAG=1)
     !!! call NC_Def_Var(ResPath,P_ResName,'SSH','m','lon','lat','time',Double_FLAG=1)
     
     call NC_Def_Var(ResPath,U_ResName,'TAUX','N/m^2','lon','lat','time',DOUBLE_FLAG=1)
     call NC_Def_Var(ResPath,V_ResName,'TAUY','N/m^2','lon','lat','time',DOUBLE_FLAG=1)
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call NC_Def_Var(ResPath,P_ResName,'Res_AbsIT','time step','time',Double_FLAG=1);
     call NC_Def_Var(ResPath,U_ResName,'Res_AbsIT','time step','time',Double_FLAG=1);
     call NC_Def_Var(ResPath,V_ResName,'Res_AbsIT','time step','time',Double_FLAG=1);
     
     call NC_Def_Var(ResPath,P_ResName,'CYear','year','time',Double_FLAG=1);
     call NC_Def_Var(ResPath,U_ResName,'CYear','year','time',Double_FLAG=1);
     call NC_Def_Var(ResPath,V_ResName,'CYear','year','time',Double_FLAG=1);
     
     call NC_Def_Var(ResPath,P_ResName,'CMonth','month','time',Double_FLAG=1);
     call NC_Def_Var(ResPath,U_ResName,'CMonth','month','time',Double_FLAG=1);
     call NC_Def_Var(ResPath,V_ResName,'CMonth','month','time',Double_FLAG=1);
     
     call NC_Def_Var(ResPath,P_ResName,'CDay','day','time',Double_FLAG=1);
     call NC_Def_Var(ResPath,U_ResName,'CDay','day','time',Double_FLAG=1);
     call NC_Def_Var(ResPath,V_ResName,'CDay','day','time',Double_FLAG=1);
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call NC_Write_One_IT(ResPath,P_ResName,'Co_Pres',OcnStat.Co_Pres,1);
     call NC_Write_One_IT(ResPath,U_ResName,'Co_Ucur',OcnStat.Co_Ucur,1);
     call NC_Write_One_IT(ResPath,V_ResName,'Co_Vcur',OcnStat.Co_Vcur,1);
     
     call NC_Write_One_IT(ResPath,P_ResName,'SST',OcnStat.SST,1);
     call NC_Write_One_IT(ResPath,P_ResName,'TSUB',OcnStat.TSUB,1);
     !!! call NC_Write_One_IT(ResPath,P_ResName,'SSH',OcnStat.SSH,1);
     
     call NC_Write_One_IT(ResPath,U_ResName,'TAUX',SrfFlux.TAUX,1);
     call NC_Write_One_IT(ResPath,V_ResName,'TAUY',SrfFlux.TAUY,1);
     
     !!! BE CAREFUL : TimeTick.Run_AbsIT >> Res_AbsIT
     call NC_Write_One_IT(ResPath,P_ResName,'Res_AbsIT',dble(TimTick.Run_AbsIT),1);
     call NC_Write_One_IT(ResPath,U_ResName,'Res_AbsIT',dble(TimTick.Run_AbsIT),1);
     call NC_Write_One_IT(ResPath,V_ResName,'Res_AbsIT',dble(TimTick.Run_AbsIT),1);
     
     !!! CDay, CMonth, CDay into ReStart File
     call NC_Write_One_IT(ResPath,P_ResName,'CYear',dble(TimTick.CYear),1);
     call NC_Write_One_IT(ResPath,U_ResName,'CYear',dble(TimTick.CYear),1);
     call NC_Write_One_IT(ResPath,V_ResName,'CYear',dble(TimTick.CYear),1);
     
     call NC_Write_One_IT(ResPath,P_ResName,'CMonth',dble(TimTick.CMonth),1);
     call NC_Write_One_IT(ResPath,U_ResName,'CMonth',dble(TimTick.CMonth),1);
     call NC_Write_One_IT(ResPath,V_ResName,'CMonth',dble(TimTick.CMonth),1);
     
     call NC_Write_One_IT(ResPath,P_ResName,'CDay',dble(TimTick.CDay),1);
     call NC_Write_One_IT(ResPath,U_ResName,'CDay',dble(TimTick.CDay),1);
     call NC_Write_One_IT(ResPath,V_ResName,'CDay',dble(TimTick.CDay),1);
     
     !!! make sure not empty file output 
     write(*,*),'SSHA_N4  = '//trim(Num2Str(OcnStat.SSHA_N4))//' SSHA_N34  = '//trim(Num2Str(OcnStat.SSHA_N34))//&
                ' SSHA_N3  = '//trim(Num2Str(OcnStat.SSHA_N3));
     write(*,*),''
     
     write(*,*),'UA_N4    = '//trim(Num2Str(OcnStat.UA_N4))//  ' UA_N34    = '//trim(Num2Str(OcnStat.UA_N34))//&
                ' UA_N3    = '//trim(Num2str(OcnStat.UA_N3));
     write(*,*),''
     
     write(*,*),'SSTA_N4  = '//trim(Num2Str(OcnStat.SSTA_N4))//' SSTA_N34  = '//trim(Num2Str(OcnStat.SSTA_N34))//&
                ' SSTA_N3  = '//trim(Num2Str(OcnStat.SSTA_N3));
     write(*,*),''
     
     
     write(*,*),'ReStart Files Writing Completed'
     write(*,*),'------------------------------------------------------------------------------'
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE Write_ReStart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Read_ReStart(ResPath,GrdInfo,OcnStat,TimTick,SrfFlux);

     TYPE (OceanStat) :: OcnStat
     TYPE (GridInfo) :: GrdInfo
     character(len=*) :: ResPath
     TYPE (TimeTick):: TimTick
     TYPE(SurfaceFlux) :: SrfFlux
     
     real*8 :: db_Res_AbsIT, db_CYear, db_CMonth, db_CDay
     !!! double precision absolute integration IT from one Run for restart
     
     character(len=300) :: P_ResName = 'ReStart_PGrid.nc';
     character(len=300) :: U_ResName = 'ReStart_UGrid.nc';
     character(len=300) :: V_ResName = 'ReStart_VGrid.nc';
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     write(*,*)'------------------------------------------------------------------------------'
     write(*,*),'Reading ReStart...'
     
     call NC_Read_One_IT(ResPath,P_ResName,'Co_Pres',OcnStat.Co_Pres,1,'N');
     call NC_Read_One_IT(ResPath,U_ResName,'Co_Ucur',OcnStat.Co_Ucur,1,'N');
     call NC_Read_One_IT(ResPath,V_ResName,'Co_Vcur',OcnStat.Co_Vcur,1,'N');
     
     call NC_Read_One_IT(ResPath,P_ResName,'SST',OcnStat.SST,1,'N');
     call NC_Read_One_IT(ResPath,P_ResName,'TSUB',OcnStat.TSUB,1,'N');
     !!! call NC_Read_One_IT(ResPath,P_ResName,'SSH',OcnStat.SSH,1,'N');
     
     call NC_Read_One_IT(ResPath,U_ResName,'TAUX',SrfFlux.TAUX,1,'N');
     call NC_Read_One_IT(ResPath,V_ResName,'TAUY',SrfFlux.TAUY,1,'N');
     

     !!! only read once from P_ResName
     call NC_Read_One_IT(ResPath,P_ResName,'Res_AbsIT',db_Res_AbsIT,1,'N');
     call NC_Read_One_IT(ResPath,P_ResName,'CYear',db_CYear,1,'N');
     call NC_Read_One_IT(ResPath,P_ResName,'CMonth',db_CMonth,1,'N');
     call NC_Read_One_IT(ResPath,P_ResName,'CDay',db_CDay,1,'N');
    
     TimTick.Res_AbsIT = int8(db_Res_AbsIT);
     TimTick.CYear     = int8(db_CYear);
     TimTick.CMonth    = int8(db_CMonth);
     TimTick.CDay      = int8(db_CDay);
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!! Clip_Run, directly assign days of run 
     IF ( Clip_Run_OnOff == 1 ) THEN
     
        TimTick.Run_AbsNT   = Clip_Run_NDay * (Day_NT) + TimTick.Res_AbsIT
        TimTick.Run_AbsNDay = Clip_Run_NDay + ceiling(dble(TimTick.Res_AbsIT)/dble(Day_NT))
    
     END IF
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
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
     
     write(*,*),'ReStart Reading Completed'
     write(*,*),'-----------------------------------------------------------------------------'
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE Read_ReStart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ReInit_Stat(GrdInfo,OcnStat,TimTick,SrfFlux);

     TYPE (OceanStat) :: OcnStat
     TYPE (GridInfo) :: GrdInfo
     TYPE (TimeTick):: TimTick
     TYPE(SurfaceFlux) :: SrfFlux
     
     OcnStat.Co_Pres = 0.0d0;
     OcnStat.Co_Ucur = 0.0d0;
     OcnStat.Co_Vcur = 0.0d0;
     
     OcnStat.SST = 0.0d0;
     OcnStat.TSUB = 0.0d0;
     
     SrfFlux.TAUX = 0.0d0;
     SrfFlux.TAUY = 0.0d0;

     
END SUBROUTINE ReInit_Stat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Create_OutPut(OutPath,GrdInfo,OcnStat,TimTick,AllFrc_OcnStat,SrfFlux);

     TYPE (OceanStat) :: OcnStat
     TYPE (GridInfo) :: GrdInfo
     character(len=*) :: OutPath
     TYPE(TimeTick) :: TimTick  
     TYPE(AllFrc_OceanStat) :: AllFrc_OcnStat
     TYPE(SurfaceFlux) :: SrfFlux
     
     character(len=300) :: P_OutName
     character(len=300) :: U_OutName 
     character(len=300) :: V_OutName 
     
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
     
     !!! to prevent mixjudge of CDO (sub half-day)
     out_time_tick = out_time_tick - 0.5d0
     
     !!! time_unit
     time_unit = &
     'days since '//trim(Math_YMD_To_Date(Start_CYear,Start_CMonth,Start_CDay))//' 00:00:0.0';
     
     !!! P_OutName, U_OutName, V_OutName
     P_OutName = 'PGrid_'//trim(TimTick.Out_Date)//'_'//trim(Num2Str(OutPut_Day_Span))//'day_'// &
                 trim(Num2Str(PutData_N_Day))//'day_R'//trim(Num2Str(TimTick.Run_Cycle_ID))//'.nc'
     
     U_OutName = 'UGrid_'//trim(TimTick.Out_Date)//'_'//trim(Num2Str(OutPut_Day_Span))//'day_'//&
                 trim(Num2Str(PutData_N_Day))//'day_R'//trim(Num2Str(TimTick.Run_Cycle_ID))//'.nc'
     
     V_OutName = 'VGrid_'//trim(TimTick.Out_Date)//'_'//trim(Num2Str(OutPut_Day_Span))//'day_'//&
                 trim(Num2Str(PutData_N_Day))//'day_R'//trim(Num2Str(TimTick.Run_Cycle_ID))//'.nc'
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     write(*,*),'----------------------------------------------------------------------------'
     write(*,*),'Creating OutPut Files ...'

     call  NC_Create_File( OutPath,P_OutName, &
                           GrdInfo.PLon, GrdInfo.PLat, GrdInfo.IMode, GrdInfo.ULDepth, out_time_tick, &
                           'lon', 'lat','imode','depth','time',time_unit, 'X', Leap_Calendar_OnOff)
    
     IF(Clip_Run_OnOff == 0) THEN
     
     call  NC_Create_File( OutPath,U_OutName, &
                           GrdInfo.ULon, GrdInfo.ULat, GrdInfo.IMode, GrdInfo.ULDepth, out_time_tick, &
                           'lon', 'lat', 'imode','depth','time',time_unit, 'X', Leap_Calendar_OnOff)
    
     call  NC_Create_File( OutPath,V_OutName, &
                           GrdInfo.VLon, GrdInfo.VLat, GrdInfo.IMode, GrdInfo.ULDepth, out_time_tick, &
                           'lon', 'lat', 'imode','depth','time',time_unit,'X', Leap_Calendar_OnOff)
     
     END IF !!!  IF(Clip_Run_OnOff == 0)
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!! OutPut timetick information
     
     IF(Clip_Run_OnOff == 0)THEN
     
     call NC_Def_Var(OutPath,P_OutName,'Res_AbsIT','time step','time',Double_FLAG=1);
     call NC_Def_Var(OutPath,U_OutName,'Res_AbsIT','time step','time',Double_FLAG=1);
     call NC_Def_Var(OutPath,V_OutName,'Res_AbsIT','time step','time',Double_FLAG=1);
     
     call NC_Def_Var(OutPath,P_OutName,'CYear','year','time',Double_FLAG=1);
     call NC_Def_Var(OutPath,U_OutName,'CYear','year','time',Double_FLAG=1);
     call NC_Def_Var(OutPath,V_OutName,'CYear','year','time',Double_FLAG=1);
     
     call NC_Def_Var(OutPath,P_OutName,'CMonth','month','time',Double_FLAG=1);
     call NC_Def_Var(OutPath,U_OutName,'CMonth','month','time',Double_FLAG=1);
     call NC_Def_Var(OutPath,V_OutName,'CMonth','month','time',Double_FLAG=1);
     
     call NC_Def_Var(OutPath,P_OutName,'CDay','day','time',Double_FLAG=1);
     call NC_Def_Var(OutPath,U_OutName,'CDay','day','time',Double_FLAG=1);
     call NC_Def_Var(OutPath,V_OutName,'CDay','day','time',Double_FLAG=1);
     
     END IF !!!  IF(Clip_Run_OnOff == 0)
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     IF(Clip_Run_OnOff == 0) THEN
     
     !!! 4D (X,Y,M,T) can be used as ReStart File
     call NC_Def_Var(OutPath,P_OutName,'Co_Pres','none','lon','lat','imode','time',Double_FLAG=1)
     call NC_Def_Var(OutPath,U_OutName,'Co_Ucur','none','lon','lat','imode','time',Double_FLAG=1)
     call NC_Def_Var(OutPath,V_OutName,'Co_Vcur','none','lon','lat','imode','time',Double_FLAG=1)
     
     !!! 4D (X,Y,Z,T)
     IF (OutPut_UpperLayer_OnOff == 1) THEN
        call NC_Def_Var(OutPath,P_OutName,'UL_DzDt','m/s','lon','lat','depth','time')
        call NC_Def_Var(OutPath,U_OutName,'UL_Ucur','m/s','lon','lat','depth','time')
        call NC_Def_Var(OutPath,V_OutName,'UL_Vcur','m/s','lon','lat','depth','time')
     END IF
     
     !!! 3D (X,Y,T)
     call NC_Def_Var(OutPath,U_OutName,'Usrf','m/s','lon','lat','time')
     call NC_Def_Var(OutPath,V_OutName,'Vsrf','m/s','lon','lat','time')
     
     call NC_Def_Var(OutPath,P_OutName,'Wsrf','m/s','lon','lat','time')
     call NC_Def_Var(OutPath,P_OutName,'SSH','m','lon','lat','time')
     
     !!! TAUX and TAUY can be used as ReStart File
     call NC_Def_Var(OutPath,U_OutName,'TAUX','N/m^2','lon','lat','time',DOUBLE_FLAG=1)
     call NC_Def_Var(OutPath,V_OutName,'TAUY','N/m^2','lon','lat','time',DOUBLE_FLAG=1)
     
     call NC_Def_Var(OutPath,U_OutName,'TAUX_Test','N/m^2','lon','lat','time')
     call NC_Def_Var(OutPath,V_OutName,'TAUY_Test','N/m^2','lon','lat','time')
     
     END IF !!! IF(Clip_Run_OnOff == 0)
     
    
     !!! SSTA tendency 3D (X,Y,T)
     IF (Solve_SSTA_OnOff == 1 ) THEN
     
         IF(Clip_Run_OnOff == 0) THEN
         
         !!! SST and TSUB can be used as ReStart File
         call NC_Def_Var(OutPath,P_OutName,'SST','dC','lon','lat','time', DOUBLE_FLAG=1)
         call NC_Def_Var(OutPath,P_OutName,'TSUB','dC','lon','lat','time',DOUBLE_FLAG=1)
         
         END IF !!! (Clip_Run_OnOff == 0)
         
         !!! heat budget
         IF( OutPut_HeatBDG_OnOff == 1 ) THEN
     
             call NC_Def_Var(OutPath,P_OutName,'UATBX','dC/s','lon','lat','time')
             call NC_Def_Var(OutPath,P_OutName,'UBTAX','dC/s','lon','lat','time')
             call NC_Def_Var(OutPath,P_OutName,'UATAX','dC/s','lon','lat','time')
     
             call NC_Def_Var(OutPath,P_OutName,'VATBY','dC/s','lon','lat','time')
             call NC_Def_Var(OutPath,P_OutName,'VBTAY','dC/s','lon','lat','time')
             call NC_Def_Var(OutPath,P_OutName,'VATAY','dC/s','lon','lat','time')
     
             call NC_Def_Var(OutPath,P_OutName,'WATBZ','dC/s','lon','lat','time')
             call NC_Def_Var(OutPath,P_OutName,'WBTAZ','dC/s','lon','lat','time')
             call NC_Def_Var(OutPath,P_OutName,'WATAZ','dC/s','lon','lat','time')
     
             call NC_Def_Var(OutPath,P_OutName,'SSTA_RESQ','dC/s','lon','lat','time')
         
         END IF !!! OutPut_HeatBDG_OnOff
         
         
         !!! heat budget
         IF( OutPut_Simple_HeatBDG_OnOff == 1 ) THEN
         
             call NC_Def_Var(OutPath,P_OutName,'ZAFK','dC/s','lon','lat','time')
             call NC_Def_Var(OutPath,P_OutName,'THFK','dC/s','lon','lat','time')
             call NC_Def_Var(OutPath,P_OutName,'RESFK','dC/s','lon','lat','time')
             call NC_Def_Var(OutPath,P_OutName,'NUDGE','dC/s','lon','lat','time')
         
         
         END IF !!! OutPut_Simple_HeatBDG_OnOff
         
         
     
     END IF  !!! Solve_SSTA_OnOff
     
     
     !!! Time series
     call NC_Def_Var(OutPath,P_OutName,'SSHA_N4','m','time')
     call NC_Def_Var(OutPath,P_OutName,'SSHA_N34','m','time')
     call NC_Def_Var(OutPath,P_OutName,'SSHA_N3','m','time')
     
     call NC_Def_Var(OutPath,P_OutName,'TSUB_N4','dC','time')
     call NC_Def_Var(OutPath,P_OutName,'TSUB_N34','dC','time')
     call NC_Def_Var(OutPath,P_OutName,'TSUB_N3','dC','time')
     
     call NC_Def_Var(OutPath,P_OutName,'SSTA_N4','dC','time')
     call NC_Def_Var(OutPath,P_OutName,'SSTA_N34','dC','time')
     call NC_Def_Var(OutPath,P_OutName,'SSTA_N3','dC','time')
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF(Clip_Run_OnOff == 0)THEN
     
     call NC_Def_Var(OutPath,U_OutName,'UA_N4','m/s','time')
     call NC_Def_Var(OutPath,U_OutName,'UA_N34','m/s','time')
     call NC_Def_Var(OutPath,U_OutName,'UA_N3','m/s','time')
     
     call NC_Def_Var(OutPath,U_OutName,'TAUX_N4','N/m^2','time')
     call NC_Def_Var(OutPath,U_OutName,'TAUX_N34','N/m^2','time')
     call NC_Def_Var(OutPath,U_OutName,'TAUX_N3','N/m^2','time')
     
     call NC_Def_Var(OutPath,U_OutName,'TAUX_N4_Test','N/m^2','time')
     call NC_Def_Var(OutPath,U_OutName,'TAUX_N34_Test','N/m^2','time')
     call NC_Def_Var(OutPath,U_OutName,'TAUX_N3_Test','N/m^2','time')
     
     ENDIF !!! IF(Clip_Run_OnOff == 0)
    
     write(*,*),'OutPut Files Creating Completed'
     write(*,*),'------------------------------------------------------------------------------'
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE Create_OutPut
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Create_MonMe_OutPut(OutPath,GrdInfo,OcnStat,TimTick,AllFrc_OcnStat,SrfFlux);

     TYPE (OceanStat) :: OcnStat
     TYPE (GridInfo) :: GrdInfo
     character(len=*) :: OutPath
     TYPE(TimeTick) :: TimTick  
     TYPE(AllFrc_OceanStat) :: AllFrc_OcnStat
     TYPE(SurfaceFlux) :: SrfFlux
     
     character(len=300) :: P_OutName
     character(len=300) :: U_OutName 
     character(len=300) :: V_OutName 
     
     integer*8 :: IT
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
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
     
     !!! P_OutName, U_OutName, V_OutName
     P_OutName = 'MonMe_PGrid_'//trim(TimTick.MonMe_Out_Date)//'_'//&
                  trim(Num2Str(MonMe_OutPut_Mon_Span))//'mon_R'//trim(Num2Str(TimTick.Run_Cycle_ID))//'.nc'
     
     U_OutName = 'MonMe_UGrid_'//trim(TimTick.MonMe_Out_Date)//'_'//&
                  trim(Num2Str(MonMe_OutPut_Mon_Span))//'mon_R'//trim(Num2Str(TimTick.Run_Cycle_ID))//'.nc'
     
     V_OutName = 'MonMe_VGrid_'//trim(TimTick.MonMe_Out_Date)//'_'//&
                  trim(Num2Str(MonMe_OutPut_Mon_Span))//'mon_R'//trim(Num2Str(TimTick.Run_Cycle_ID))//'.nc'
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     write(*,*),'----------------------------------------------------------------------------'
     write(*,*),'Creating Monthly Mean OutPut Files ...'

     call  NC_Create_File( OutPath,P_OutName, &
                           GrdInfo.PLon, GrdInfo.PLat, GrdInfo.IMode, GrdInfo.ULDepth, out_time_tick, &
                           'lon', 'lat','imode','depth','time',time_unit, 'X', Leap_Calendar_OnOff)
    
     call  NC_Create_File( OutPath,U_OutName, &
                           GrdInfo.ULon, GrdInfo.ULat, GrdInfo.IMode, GrdInfo.ULDepth, out_time_tick, &
                           'lon', 'lat', 'imode','depth','time',time_unit, 'X', Leap_Calendar_OnOff)
    
     call  NC_Create_File( OutPath,V_OutName, &
                           GrdInfo.VLon, GrdInfo.VLat, GrdInfo.IMode, GrdInfo.ULDepth, out_time_tick, &
                           'lon', 'lat', 'imode','depth','time',time_unit,'X', Leap_Calendar_OnOff)
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
     !!! 3D (X,Y,T)
     call NC_Def_Var(OutPath,U_OutName,'Usrf','m/s','lon','lat','time')
     call NC_Def_Var(OutPath,V_OutName,'Vsrf','m/s','lon','lat','time')
     
     call NC_Def_Var(OutPath,P_OutName,'Wsrf','m/s','lon','lat','time')
     call NC_Def_Var(OutPath,P_OutName,'SSH','m','lon','lat','time')
     
     !!! TAUX and TAUY can be used as ReStart File
     call NC_Def_Var(OutPath,U_OutName,'TAUX','N/m^2','lon','lat','time')
     call NC_Def_Var(OutPath,V_OutName,'TAUY','N/m^2','lon','lat','time')
     
     call NC_Def_Var(OutPath,U_OutName,'TAUX_Test','N/m^2','lon','lat','time')
     call NC_Def_Var(OutPath,V_OutName,'TAUY_Test','N/m^2','lon','lat','time')
     

     !!! SSTA tendency 3D (X,Y,T)
     IF (Solve_SSTA_OnOff == 1 ) THEN
     
         !!! SST and TSUB can be used as ReStart File
         call NC_Def_Var(OutPath,P_OutName,'SST','dC','lon','lat','time')
         call NC_Def_Var(OutPath,P_OutName,'TSUB','dC','lon','lat','time')
         
         !!! heat budget
         !!! IF( OutPut_Simple_HeatBDG_OnOff == 1 ) THEN
         
             call NC_Def_Var(OutPath,P_OutName,'ZAFK','dC/s','lon','lat','time')
             call NC_Def_Var(OutPath,P_OutName,'THFK','dC/s','lon','lat','time')
             call NC_Def_Var(OutPath,P_OutName,'RESFK','dC/s','lon','lat','time')
             call NC_Def_Var(OutPath,P_OutName,'NUDGE','dC/s','lon','lat','time')
         
         
         !!! END IF !!! OutPut_Simple_HeatBDG_OnOff
     
     END IF  !!! Solve_SSTA_OnOff
     
     call NC_Def_Var(OutPath,P_OutName,'SSTA_N4','dC','time')
     call NC_Def_Var(OutPath,P_OutName,'SSTA_N34','dC','time')
     call NC_Def_Var(OutPath,P_OutName,'SSTA_N3','dC','time')
    
     write(*,*),'Monthly Mean OutPut Files Creating Completed'
     write(*,*),'------------------------------------------------------------------------------'
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
END SUBROUTINE Create_MonMe_OutPut
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Write_OutPut(OutPath,GrdInfo,OcnStat,TimTick,AllFrc_OcnStat,SrfFlux,WriteIT);

     TYPE (OceanStat) :: OcnStat
     TYPE (GridInfo) :: GrdInfo
     character(len=*) :: OutPath
     TYPE(TimeTick) :: TimTick 
     TYPE(AllFrc_OceanStat) :: AllFrc_OcnStat
     TYPE(SurfaceFlux) :: SrfFlux
     
     integer*4 :: WriteIT !!!! controlled by Ctrl_OutPut_Stream
     
     character(len=300) :: P_OutName
     character(len=300) :: U_OutName 
     character(len=300) :: V_OutName 
     
     !!! P_OutName, U_OutName, V_OutName
     P_OutName = 'PGrid_'//trim(TimTick.Out_Date)//'_'//trim(Num2Str(OutPut_Day_Span))//'day_'// &
                 trim(Num2Str(PutData_N_Day))//'day_R'//trim(Num2Str(TimTick.Run_Cycle_ID))//'.nc'
     
     U_OutName = 'UGrid_'//trim(TimTick.Out_Date)//'_'//trim(Num2Str(OutPut_Day_Span))//'day_'//&
                 trim(Num2Str(PutData_N_Day))//'day_R'//trim(Num2Str(TimTick.Run_Cycle_ID))//'.nc'
     
     V_OutName = 'VGrid_'//trim(TimTick.Out_Date)//'_'//trim(Num2Str(OutPut_Day_Span))//'day_'//&
                 trim(Num2Str(PutData_N_Day))//'day_R'//trim(Num2Str(TimTick.Run_Cycle_ID))//'.nc'
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
     IF(Clip_Run_OnOff == 0)THEN
      
     !!! BE CAREFUL : TimeTick.Run_AbsIT >> Res_AbsIT
     call NC_Write_One_IT(OutPath,P_OutName,'Res_AbsIT',dble(TimTick.Run_AbsIT),WriteIT);
     call NC_Write_One_IT(OutPath,U_OutName,'Res_AbsIT',dble(TimTick.Run_AbsIT),WriteIT);
     call NC_Write_One_IT(OutPath,V_OutName,'Res_AbsIT',dble(TimTick.Run_AbsIT),WriteIT);
     
     !!! CDay, CMonth, CDay into OutPut File
     call NC_Write_One_IT(OutPath,P_OutName,'CYear',dble(TimTick.CYear),WriteIT);
     call NC_Write_One_IT(OutPath,U_OutName,'CYear',dble(TimTick.CYear),WriteIT);
     call NC_Write_One_IT(OutPath,V_OutName,'CYear',dble(TimTick.CYear),WriteIT);
     
     call NC_Write_One_IT(OutPath,P_OutName,'CMonth',dble(TimTick.CMonth),WriteIT);
     call NC_Write_One_IT(OutPath,U_OutName,'CMonth',dble(TimTick.CMonth),WriteIT);
     call NC_Write_One_IT(OutPath,V_OutName,'CMonth',dble(TimTick.CMonth),WriteIT);
     
     call NC_Write_One_IT(OutPath,P_OutName,'CDay',dble(TimTick.CDay),WriteIT);
     call NC_Write_One_IT(OutPath,U_OutName,'CDay',dble(TimTick.CDay),WriteIT);
     call NC_Write_One_IT(OutPath,V_OutName,'CDay',dble(TimTick.CDay),WriteIT);
     
     END IF !!!  IF(Clip_Run_OnOff == 0)
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     
     IF(Clip_Run_OnOff == 0)THEN
      
     !!! coefficient of baroclinic mode
     call NC_Write_One_IT(OutPath,P_OutName,'Co_Pres',OcnStat.Co_Pres,WriteIT);
     call NC_Write_One_IT(OutPath,U_OutName,'Co_Ucur',OcnStat.Co_Ucur,WriteIT);
     call NC_Write_One_IT(OutPath,V_OutName,'Co_Vcur',OcnStat.Co_Vcur,WriteIT);
     
     !!! surface field
     call NC_Write_One_IT(OutPath,U_OutName,'Usrf',OcnStat.Usrf,WriteIT);
     call NC_Write_One_IT(OutPath,V_OutName,'Vsrf',OcnStat.Vsrf,WriteIT);
     
     call NC_Write_One_IT(OutPath,P_OutName,'Wsrf',OcnStat.Wsrf,WriteIT);
     call NC_Write_One_IT(OutPath,P_OutName,'SSH',OcnStat.SSH,WriteIT);
     
     call NC_Write_One_IT(OutPath,U_OutName,'TAUX',SrfFlux.TAUX,WriteIT);
     call NC_Write_One_IT(OutPath,U_OutName,'TAUX_Test',SrfFlux.TAUX_Test,WriteIT);
     
     call NC_Write_One_IT(OutPath,V_OutName,'TAUY',SrfFlux.TAUY,WriteIT);
     call NC_Write_One_IT(OutPath,V_OutName,'TAUY_Test',SrfFlux.TAUY_Test,WriteIT);
     
     END IF !!! IF(Clip_Run_OnOff == 0)
     
     
     !!! SSTA tendency
     IF (Solve_SSTA_OnOff == 1 ) THEN
     
         IF(Clip_Run_OnOff == 0) THEN
     
         call NC_Write_One_IT(OutPath,P_OutName,'SST',OcnStat.SST,WriteIT);
         call NC_Write_One_IT(OutPath,P_OutName,'TSUB',OcnStat.TSUB,WriteIT);
         
         END IF !!! IF(Clip_Run_OnOff == 0) THEN
         
         !!! Heat budget
         IF( OutPut_HeatBDG_OnOff == 1 ) THEN
     
             call NC_Write_One_IT(OutPath,P_OutName,'UATBX',AllFrc_OcnStat.UATBX,WriteIT);
             call NC_Write_One_IT(OutPath,P_OutName,'UBTAX',AllFrc_OcnStat.UBTAX,WriteIT);
             call NC_Write_One_IT(OutPath,P_OutName,'UATAX',AllFrc_OcnStat.UATAX,WriteIT);
     
             call NC_Write_One_IT(OutPath,P_OutName,'VATBY',AllFrc_OcnStat.VATBY,WriteIT);
             call NC_Write_One_IT(OutPath,P_OutName,'VBTAY',AllFrc_OcnStat.VBTAY,WriteIT);
             call NC_Write_One_IT(OutPath,P_OutName,'VATAY',AllFrc_OcnStat.VATAY,WriteIT);
     
             call NC_Write_One_IT(OutPath,P_OutName,'WATBZ',AllFrc_OcnStat.WATBZ,WriteIT);
             call NC_Write_One_IT(OutPath,P_OutName,'WBTAZ',AllFrc_OcnStat.WBTAZ,WriteIT);
             call NC_Write_One_IT(OutPath,P_OutName,'WATAZ',AllFrc_OcnStat.WATAZ,WriteIT);
     
             call NC_Write_One_IT(OutPath,P_OutName,'SSTA_RESQ',AllFrc_OcnStat.SSTA_RESQ,WriteIT);
         
         END IF !!! OutPut_Simple_HeatBDG_OnOff
         
         
         !!! Simple Heat Budget
         IF( OutPut_Simple_HeatBDG_OnOff == 1 ) THEN
         
             call NC_Write_One_IT(OutPath,P_OutName,'ZAFK',AllFrc_OcnStat.ZAFK,WriteIT);
             call NC_Write_One_IT(OutPath,P_OutName,'THFK',AllFrc_OcnStat.THFK,WriteIT);
             call NC_Write_One_IT(OutPath,P_OutName,'RESFK',AllFrc_OcnStat.RESFK,WriteIT);
             call NC_Write_One_IT(OutPath,P_OutName,'NUDGE',AllFrc_OcnStat.NUDGE,WriteIT);
         
         END IF !!! OutPut_Simple_HeatBDG_OnOff
         
     
     END IF !!! Solve_SSTA_OnOff
     
     
     !!! time series
     call NC_Write_One_IT(OutPath,P_OutName,'SSHA_N4',OcnStat.SSHA_N4,WriteIT);
     call NC_Write_One_IT(OutPath,P_OutName,'SSHA_N34',OcnStat.SSHA_N34,WriteIT);
     call NC_Write_One_IT(OutPath,P_OutName,'SSHA_N3',OcnStat.SSHA_N3,WriteIT);
     
     call NC_Write_One_IT(OutPath,P_OutName,'TSUB_N4',OcnStat.TSUB_N4,WriteIT);
     call NC_Write_One_IT(OutPath,P_OutName,'TSUB_N34',OcnStat.TSUB_N34,WriteIT);
     call NC_Write_One_IT(OutPath,P_OutName,'TSUB_N3',OcnStat.TSUB_N3,WriteIT);
     
     call NC_Write_One_IT(OutPath,P_OutName,'SSTA_N4',OcnStat.SSTA_N4,WriteIT);
     call NC_Write_One_IT(OutPath,P_OutName,'SSTA_N34',OcnStat.SSTA_N34,WriteIT);
     call NC_Write_One_IT(OutPath,P_OutName,'SSTA_N3',OcnStat.SSTA_N3,WriteIT);
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF(Clip_Run_OnOff == 0)THEN
     
     call NC_Write_One_IT(OutPath,U_OutName,'UA_N4',OcnStat.UA_N4,WriteIT);
     call NC_Write_One_IT(OutPath,U_OutName,'UA_N34',OcnStat.UA_N34,WriteIT);
     call NC_Write_One_IT(OutPath,U_OutName,'UA_N3',OcnStat.UA_N3,WriteIT);
     
     call NC_Write_One_IT(OutPath,U_OutName,'TAUX_N4',SrfFlux.TAUX_N4,WriteIT);
     call NC_Write_One_IT(OutPath,U_OutName,'TAUX_N34',SrfFlux.TAUX_N34,WriteIT);
     call NC_Write_One_IT(OutPath,U_OutName,'TAUX_N3',SrfFlux.TAUX_N3,WriteIT);
     
     call NC_Write_One_IT(OutPath,U_OutName,'TAUX_N4_Test',SrfFlux.TAUX_N4_Test,WriteIT);
     call NC_Write_One_IT(OutPath,U_OutName,'TAUX_N34_Test',SrfFlux.TAUX_N34_Test,WriteIT);
     call NC_Write_One_IT(OutPath,U_OutName,'TAUX_N3_Test',SrfFlux.TAUX_N3_Test,WriteIT);
     
     ENDIF !!! IF(Clip_Run_OnOff == 0)
     
     !!! upperlayer field
     IF (OutPut_UpperLayer_OnOff == 1) THEN
     
        call NC_Write_One_IT(OutPath,P_OutName,'UL_DzDt',OcnStat.UL_DzDt,WriteIT);
        call NC_Write_One_IT(OutPath,U_OutName,'UL_Ucur',OcnStat.UL_Ucur,WriteIT);
        call NC_Write_One_IT(OutPath,V_OutName,'UL_Vcur',OcnStat.UL_Vcur,WriteIT);
     
     END IF
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF ( Write_OutPut_Check == 1 ) THEN
     
     write(*,*),'------------------------------------------------------------------------------'
     write(*,*),'Writing OutPut Files ...'
     
     write(6,101),WriteIT
     101 format(' WriteIT = ',I8);  
     write(*,*),'------------------------------------------------------------------------------'
     
     !!! Make sure not empty file output & check if openmp works correctly
     !!! surface flux series
     !!!write(*,*),'TAUX_N4       = '//trim(Num2Str(SrfFlux.TAUX_N4))//&
     !!!           '  TAUX_N34       = '//trim(Num2Str(SrfFlux.TAUX_N34))//&
     !!!           '  TAUX_N3        = '//trim(Num2Str(SrfFlux.TAUX_N3));
     !!!write(*,*),''
     !!!
     !!!write(*,*),'TAUX_N4_Test  = '//trim(Num2Str(SrfFlux.TAUX_N4_Test))//&
     !!!           '  TAUX_N34_Test  = '//trim(Num2Str(SrfFlux.TAUX_N34_Test))//&
     !!!           '  TAUX_N3_Test   = '//trim(Num2Str(SrfFlux.TAUX_N3_Test));
     !!!write(*,*),''
     
     !!!! ocean stat series
     write(*,*),'SSHA_N4       = '//trim(Num2Str(OcnStat.SSHA_N4))//&
                '  SSHA_N34       = '//trim(Num2Str(OcnStat.SSHA_N34))//&
                '  SSHA_N3        = '//trim(Num2Str(OcnStat.SSHA_N3));
     write(*,*),''
     
     !!!write(*,*),'UA_N4         = '//trim(Num2Str(OcnStat.UA_N4))//&
     !!!           '  UA_N34         = '//trim(Num2Str(OcnStat.UA_N34))//&
     !!!           '  UA_N3          = '//trim(Num2str(OcnStat.UA_N3));
     !!!write(*,*),''
     !!!
     !!!write(*,*),'TSUB_N4       = '//trim(Num2Str(OcnStat.TSUB_N4))//&
     !!!           '  TSUB_N34       = '//trim(Num2Str(OcnStat.TSUB_N34))//&
     !!!           '  TSUB_N3        = '//trim(Num2Str(OcnStat.TSUB_N3));
     !!!write(*,*),''
     
     write(*,*),'SSTA_N4       = '//trim(Num2Str(OcnStat.SSTA_N4))//&
                '  SSTA_N34       = '//trim(Num2Str(OcnStat.SSTA_N34))//&
                '  SSTA_N3        = '//trim(Num2Str(OcnStat.SSTA_N3));
     write(*,*),''
     
     !!!write(*,*),'Eq_SSTA_N4    = '//trim(Num2Str(OcnStat.Eq_SSTA_N4))//&
     !!!           '  Eq_SSTA_N34    = '//trim(Num2Str(OcnStat.Eq_SSTA_N34))//&
     !!!           '  Eq_SSTA_N3     = '//trim(Num2Str(OcnStat.Eq_SSTA_N3));
     !!!write(*,*),''
     
     
     write(*,*),'OutPut Files Writing Completed'
     write(*,*),'------------------------------------------------------------------------------'
     
     END IF !!! Write_OutPut_Check
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE Write_OutPut
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Write_MonMe_OutPut(OutPath,GrdInfo,OcnStat,TimTick,AllFrc_OcnStat,SrfFlux,WriteIT);

     TYPE (OceanStat) :: OcnStat
     TYPE (GridInfo) :: GrdInfo
     character(len=*) :: OutPath
     TYPE(TimeTick) :: TimTick 
     TYPE(AllFrc_OceanStat) :: AllFrc_OcnStat
     TYPE(SurfaceFlux) :: SrfFlux
     
     integer*4 :: WriteIT !!!! controlled by Ctrl_OutPut_Stream
     
     character(len=300) :: P_OutName
     character(len=300) :: U_OutName 
     character(len=300) :: V_OutName 
     
     !!! P_OutName, U_OutName, V_OutName
     P_OutName = 'MonMe_PGrid_'//trim(TimTick.MonMe_Out_Date)//'_'//&
                 trim(Num2Str(MonMe_OutPut_Mon_Span))//'mon_R'//trim(Num2Str(TimTick.Run_Cycle_ID))//'.nc'
     
     U_OutName = 'MonMe_UGrid_'//trim(TimTick.MonMe_Out_Date)//'_'//&
                 trim(Num2Str(MonMe_OutPut_Mon_Span))//'mon_R'//trim(Num2Str(TimTick.Run_Cycle_ID))//'.nc'
     
     V_OutName = 'MonMe_VGrid_'//trim(TimTick.MonMe_Out_Date)//'_'//&
                 trim(Num2Str(MonMe_OutPut_Mon_Span))//'mon_R'//trim(Num2Str(TimTick.Run_Cycle_ID))//'.nc'
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!! surface field
     call NC_Write_One_IT(OutPath,U_OutName,'Usrf',OcnStat.MonMe_Usrf,WriteIT);
     call NC_Write_One_IT(OutPath,V_OutName,'Vsrf',OcnStat.MonMe_Vsrf,WriteIT);
     
     call NC_Write_One_IT(OutPath,P_OutName,'Wsrf',OcnStat.MonMe_Wsrf,WriteIT);
     call NC_Write_One_IT(OutPath,P_OutName,'SSH',OcnStat.MonMe_SSH,WriteIT);
     
     call NC_Write_One_IT(OutPath,U_OutName,'TAUX',SrfFlux.MonMe_TAUX,WriteIT);
     call NC_Write_One_IT(OutPath,U_OutName,'TAUX_Test',SrfFlux.MonMe_TAUX_Test,WriteIT);
     
     call NC_Write_One_IT(OutPath,V_OutName,'TAUY',SrfFlux.MonMe_TAUY,WriteIT);
     call NC_Write_One_IT(OutPath,V_OutName,'TAUY_Test',SrfFlux.MonMe_TAUY_Test,WriteIT);
     
     
     !!! SSTA tendency
     IF (Solve_SSTA_OnOff == 1 ) THEN
     
         call NC_Write_One_IT(OutPath,P_OutName,'SST',OcnStat.MonMe_SST,WriteIT);
         call NC_Write_One_IT(OutPath,P_OutName,'TSUB',OcnStat.MonMe_TSUB,WriteIT);
         
         !!! Simple Heat Budget
         !!! IF( OutPut_Simple_HeatBDG_OnOff == 1 ) THEN
         
             call NC_Write_One_IT(OutPath,P_OutName,'ZAFK',AllFrc_OcnStat.MonMe_ZAFK,WriteIT);
             call NC_Write_One_IT(OutPath,P_OutName,'THFK',AllFrc_OcnStat.MonMe_THFK,WriteIT);
             call NC_Write_One_IT(OutPath,P_OutName,'RESFK',AllFrc_OcnStat.MonMe_RESFK,WriteIT);
             call NC_Write_One_IT(OutPath,P_OutName,'NUDGE',AllFrc_OcnStat.MonMe_NUDGE,WriteIT);
         
         !!!! END IF !!! OutPut_Simple_HeatBDG_OnOff
         
     
     END IF !!! Solve_SSTA_OnOff
     
     call NC_Write_One_IT(OutPath,P_OutName,'SSTA_N4',OcnStat.MonMe_SSTA_N4,WriteIT);
     call NC_Write_One_IT(OutPath,P_OutName,'SSTA_N34',OcnStat.MonMe_SSTA_N34,WriteIT);
     call NC_Write_One_IT(OutPath,P_OutName,'SSTA_N3',OcnStat.MonMe_SSTA_N3,WriteIT);

     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF ( Write_OutPut_Check == 1 ) THEN
     
     write(*,*),'---------------------------------------------------------------------------------'
     write(*,*),'Writing Monthly Mean OutPut Files ...'
     
     write(6,101),WriteIT
     101 format(' WriteIT = ',I8);  
     write(*,*),'---------------------------------------------------------------------------------'
    
     write(*,*),'Monthly Mean OutPut Files Writing Completed'
     write(*,*),'---------------------------------------------------------------------------------'
     
     END IF !!! Write_OutPut_Check
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE Write_MonMe_OutPut
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Ctrl_ReStart_Stream(TimTick,GrdInfo, OcnStat,AllFrc_OcnStat, SrfFlux)

    TYPE(TimeTick) :: TimTick
    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat
    TYPE(AllFrc_OceanStat) :: AllFrc_OcnStat
    TYPE(SurfaceFlux) :: SrfFlux
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
    IF ( (TimTick.Res_FLAG_Count >= 1).AND.(mod(TimTick.Res_FLAG_Count, ReStart_N_Unit) == 0) )then
        
        TimTick.Res_FLAG_Count = 0 !!!! very important, be careful
    
        !!! update ReStart
        call Write_ReStart(ReStart_Path,GrdInfo,OcnStat,TimTick, SrfFlux);
          
        IF (TimTick.RES_FLAG_Total_Count == ReStart_N_Unit * ReStart_Cycles) then
        
            IF ( TimTick.Run_AbsIT < TimTick.Run_AbsNT ) then 
            !!! to avoid STOP caused by ReStart at the end of one run
        
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
            
            END IF !!! TimTick.Run_AbsIT < TimTick.Run_AbsNT
            
        END IF !!! TimTick.RES_FLAG_Total_Count == ReStart_N_Unit * ReStart_Cycles
    
    END IF !!! (TimTick.Res_FLAG_Count >= 1).AND.(mod(TimTick.Res_FLAG_Count, ReStart_N_Unit) == 0)

END SUBROUTINE Ctrl_ReStart_Stream
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Ctrl_OutPut_Stream(TimTick,GrdInfo, OcnStat,AllFrc_OcnStat,SrfFlux) 
!!! should be used after Ctrl_ReStart_Stream

    TYPE(TimeTick) :: TimTick
    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat 
    TYPE(AllFrc_OceanStat) :: AllFrc_OcnStat
    TYPE(SurfaceFlux) :: SrfFlux
    
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
    
            call Create_OutPut(OutPut_Path,GrdInfo,OcnStat,TimTick,AllFrc_OcnStat, SrfFlux);
        
        ENDIF
        
    ELSEIF( Clip_Run_OnOff == 1 ) THEN
    
        IF ( TimTick.Run_AbsIT == (TimTick.Res_AbsIT + 1) ) THEN
        
            !!! update TimTick.Out_Date, help target the file to write
            TimTick.Out_Date = TimTick.Date
    
            call Create_OutPut(OutPut_Path,GrdInfo,OcnStat,TimTick,AllFrc_OcnStat, SrfFlux);
        
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
    
             call Write_OutPut(OutPut_Path,GrdInfo,OcnStat,TimTick,AllFrc_OcnStat, SrfFlux, WriteIT);
    
        END IF
    ELSEIF( Clip_Run_OnOff == 1 ) THEN
    
        IF ( (TimTick.Run_AbsIT>=1) .AND. (mod(TimTick.Run_AbsIT - TimTick.Res_AbsIT,(Day_NT)*(PutData_N_Day)) == 0) )THEN
        
            !!! be careful 
            WriteIT = int4 ( ceiling(dble(TimTick.Run_AbsIT - TimTick.Res_AbsIT)/dble(Day_NT*PutData_N_Day)));
            
            !!! write(*,*),'TimTick.Run_AbsIT=',TimTick.Run_AbsIT
        
            call Write_OutPut(OutPut_Path,GrdInfo,OcnStat,TimTick,AllFrc_OcnStat, SrfFlux, WriteIT);
    
        END IF
        
    END IF 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END IF !!!  IF ( Inst_OutPut_OnOff == 1 ) THEN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

    !!! Monthly Mean OutPut
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( MonMe_OutPut_OnOff == 1 ) THEN 
    
        !!! Create 
        IF( (TimTick.Arrive_CMonth_FLAG == 1).AND.( mod(TimTick.Run_AbsIMon, MonMe_OutPut_Mon_Span)==1 ) )THEN
    
            TimTick.MonMe_Out_Date = TimTick.Date
            
            call Create_MonMe_OutPut(OutPut_Path,GrdInfo,OcnStat,TimTick,AllFrc_OcnStat, SrfFlux);
        
        END IF
        
        !!! Write
        IF( TimTick.Reach_CMonth_FLAG == 1 )THEN
        
            WriteIT = int4( mod(TimTick.Run_AbsIMon-1, MonMe_OutPut_Mon_Span) +1 );
            
        
            call Write_MonMe_OutPut(OutPut_Path,GrdInfo,OcnStat,TimTick,AllFrc_OcnStat, SrfFlux, WriteIT);
        
        
        END IF 
    
    END IF !!! MonMe_OutPut_OnOff
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
END SUBROUTINE Ctrl_OutPut_Stream
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Ctrl_Initialization (GrdInfo,Clim_OcnStat, TimTick, &
                                OcnStat, Tmp_OcnStat, AllFrc_OcnStat, SrfFlux, SrfFlux_Hist)

    TYPE(GridInfo) :: GrdInfo
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    
    TYPE(OceanStat) :: OcnStat, Tmp_OcnStat 
    !!!! Tmp_OcnStat for matsuno iteration scheme
   
    TYPE(AllFrc_OceanStat) :: AllFrc_OcnStat
    
    TYPE(SurfaceFlux) :: SrfFlux
    TYPE(SurfaceFlux_Hist) :: SrfFlux_Hist
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF( TimTick.Run_Cycle_ID == 1 ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*),'---------------------------------------------------------------------------'
    write(*,*),'Begin Initialization, Allocating Computation Space ...'
    
    !!! read grid & clim_ocnstat
    call Read_Grid_ClimStat (Grid_Path, ClimOcn_Path, Forcing_Path, GrdInfo,Clim_OcnStat)
    
    !!! allocate compuation space
    call Alloc_OceanStat (GrdInfo, OcnStat);
    call Alloc_OceanStat (GrdInfo, Tmp_OcnStat);
    call Alloc_AllFrc_OceanStat(GrdInfo, AllFrc_OcnStat);
    
    !!! allocate sea surface flux & Hist
    call Alloc_SurfaceFlux( GrdInfo, SrfFlux)
    
    call Alloc_SurfaceFlux_Hist( GrdInfo, SrfFlux_Hist, HistNT = 366) !!! consider leap year 
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Show Key Parameters on Screen
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*),'-------------------------------------------------------------------------------' 
    write(*,*),'Run_OnOff            = '//trim(Num2Str(Run_OnOff))
    write(*,*),'ReRun_OnOff          = '//trim(Num2Str(ReRun_OnOff))
    write(*,*),'Solver_OnOff         = '//trim(Num2Str(Solver_OnOff))
    write(*,*),'OutPut_OnOff         = '//trim(Num2Str(OutPut_OnOff))
    write(*,*),'ReStart_OnOff        = '//trim(Num2Str(ReStart_OnOff))
    write(*,*),'-------------------------------------------------------------------------------' 
    write(*,*),'Solve_SSTA_OnOff     = '//trim(Num2Str(Solve_SSTA_OnOff))
    write(*,*),'Solve_TSUB_AI_OnOff  = '//trim(Num2Str(Solve_TSUB_AI_OnOff))
    write(*,*),'Solve_TAUXY_AI_OnOff = '//trim(Num2Str(Solve_TAUXY_AI_OnOff))
    write(*,*),'Complex_HeatBDG_Out  = '//trim(Num2Str(OutPut_HeatBDG_OnOff))
    write(*,*),'Simple_HeatBDG_Out   = '//trim(Num2Str(OutPut_Simple_HeatBDG_OnOff))
    write(*,*),'SSTA_Nudging_OnOff   = '//trim(Num2Str(SSTA_Nudging_OnOff))
    write(*,*),'Wind_Nudging_OnOff   = '//trim(Num2Str(Wind_Nudging_OnOff))
    write(*,*),'Eq_Direct_Assign     = '//trim(Num2Str(Eq_Direct_Assign))
    write(*,*),'-------------------------------------------------------------------------------' 
    write(*,*),'Tsub_Couple_N_Day    = '//trim(Num2Str(Tsub_Couple_N_Day))
    write(*,*),'Tau_Couple_N_Day     = '//trim(Num2Str(Tau_Couple_N_Day))
    write(*,*),'Couple_Strength      = '//trim(Num2Str(Couple_Strength))
    write(*,*),'Couple_TAUX_Bias     = '//trim(Num2Str(Couple_TAUX_Bias))
    write(*,*),'Wstr_Proj_Amp        = '//trim(Num2Str(Wstr_Proj_Amp))
    write(*,*),'TSUB_InPut_Amp       = '//trim(Num2Str(TSUB_InPut_Amp))
    write(*,*),'TSUB_OutPut_Amp      = '//trim(Num2Str(TSUB_OutPut_Amp))
    write(*,*),'PhaseLocking_OnOff   = '//trim(Num2Str(PhaseLocking_OnOff))
    write(*,*),'Phase_Difference     = '//trim(Num2Str(Phase_Difference))
    write(*,*),'-------------------------------------------------------------------------------' 
    
    write(*,*),'WATBZ_Amp   = '//trim(Num2Str(WATBZ_Amp))
    write(*,*),'WBTAZ_Amp   = '//trim(Num2Str(WBTAZ_Amp))
    write(*,*),'WATAZ_Amp   = '//trim(Num2Str(WATAZ_Amp))
    
    write(*,*),'DiffX_SSTA  = '//trim(Num2Str(DiffX_SSTA))//' (m^2/s)'
    write(*,*),'DiffY_SSTA  = '//trim(Num2Str(DiffY_SSTA))//' (m^2/s)'
    write(*,*),'DiffZ_SSTA  = '//trim(Num2Str(DiffZ_SSTA))//' (m^2/s)'
    
    write(*,*),'SSTA Damping Scale = '//trim(Num2Str( (1.0d0/Damp_SSTA)/Day_DT ) )//' (day)'
    write(*,*),'SSTA_Nudging_Scale = '//trim(Num2Str( SSTA_Nudging_Scale ))//' (day)'
    write(*,*),'Wind_Nudging_Ratio = '//trim(Num2Str( Wind_Nudging_Ratio ))
    
    write(*,*),'-------------------------------------------------------------------------------'
    write(*,*),'Initialization Completed, Computation Space Allocated ...'
    write(*,*),'-------------------------------------------------------------------------------' 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
END IF !!! TimTick.Run_Cycle_ID == 1

    !!! TimTick including FTA initialization
    call Alloc_TimeTick(TimTick)

END SUBROUTINE Ctrl_Initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Ctrl_Add_External_Forcing (GrdInfo, Clim_OcnStat, OcnStat, TimTick, SrfFlux, SrfFlux_Hist )

    TYPE(GridInfo) :: GrdInfo
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(OceanStat) :: OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(SurfaceFlux) :: SrfFlux
    TYPE(SurfaceFlux_Hist) :: SrfFlux_Hist
    integer*8 :: NX,NY
    
    integer*4 :: StartIDay, CountIDay
    integer*4 :: StartIMon, CountIMon
    
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( Forcing_DT_Unit == 'C' )THEN

       !!! YYYY.01.01 first time step, read data
       !!! IF( (TimTick.CMonth == 1).AND.(TimTick.CDay == 1).AND.(mod(TimTick.Run_AbsIT,Day_NT)==1) ) THEN
       !!!
       !!! CountIMon = int4(12);
       !!! StartIMon = int4(12*(TimTick.CYear - Start_CYear) + 1);
       !!!
       !!!
       !!! call NC_Read_Mul_IT(Forcing_Path,SSH_Test_FName,'zos',&
       !!!                    SrfFlux_Hist.SSH_Hist, StartIMon, CountIMon,'N');
       !!! 
       !!! END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSEIF ( Forcing_DT_Unit == 'D' )THEN
    
        !!! YYYY.01.01 + first time step, read external forcing data
        IF( (TimTick.CMonth == 1).AND.(TimTick.CDay == 1).AND.(mod(TimTick.Run_AbsIT,Day_NT)==1) ) THEN
    
        CountIDay = int4( 365 + (Leap_Calendar_Onoff) * Math_IsLeapYear(TimTick.CYear) );
        StartIDay = int4( TimTick.Run_AbsIDay ); 
        
        !!! pseudostress forcing
        !!! call NC_Read_Mul_IT(Forcing_Path,PsTauX_Forcing_FName,'PsTauX',&
        !!!                   SrfFlux_Hist.TAUX_Hist, StartIDay, CountIDay,'N');
        !!!
        !!! call NC_Read_Mul_IT(Forcing_Path,PsTauY_Forcing_FName,'PsTauY',&
        !!!                   SrfFlux_Hist.TAUY_Hist, StartIDay, CountIDay,'N');
        !!!
        !!!!!! PsTAUXY ==> TAUXY
        !!! SrfFlux_Hist.TAUX_Hist = CoeffDrag * RhoAtm * SrfFlux_Hist.TAUX_Hist
        !!! SrfFlux_Hist.TAUY_Hist = CoeffDrag * RhoAtm * SrfFlux_Hist.TAUY_Hist
        
        call NC_Read_Mul_IT(Forcing_Path,TauX_Forcing_FName,'TauX',&
                           SrfFlux_Hist.TAUX_Hist, StartIDay, CountIDay,'N');
        
        call NC_Read_Mul_IT(Forcing_Path,TauY_Forcing_FName,'TauY',&
                           SrfFlux_Hist.TAUY_Hist, StartIDay, CountIDay,'N');
    
        END IF
        
    ELSEIF ( Forcing_DT_Unit == 'M' )THEN
    
        !!! YYYY.01.01 first time step, read data
        IF( (TimTick.CMonth == 1).AND.(TimTick.CDay == 1).AND.(mod(TimTick.Run_AbsIT,Day_NT)==1) ) THEN
    
        CountIMon = int4(12);
        StartIMon = int4(12*(TimTick.CYear - Start_CYear) + 1); 
    
        !!!call NC_Read_Mul_IT(Forcing_Path,CMIP_TauX_Forcing_FName,'CMIP_TauX',&
        !!!                   SrfFlux_Hist.TAUX_Hist, StartIMon, CountIMon,'N');
        !!!
        !!!call NC_Read_Mul_IT(Forcing_Path,CMIP_TauY_Forcing_FName,'CMIP_TauY',&
        !!!                   SrfFlux_Hist.TAUY_Hist, StartIMon, CountIMon,'N');

        call NC_Read_Mul_IT(Forcing_Path,TauX_Forcing_FName,'TauX',&
                           SrfFlux_Hist.TAUX_Hist, StartIMon, CountIMon,'N');
        
        call NC_Read_Mul_IT(Forcing_Path,TauY_Forcing_FName,'TauY',&
                           SrfFlux_Hist.TAUY_Hist, StartIMon, CountIMon,'N');
    

        !!! PsTAUXY ==> TAUXY
        !!! SrfFlux_Hist.TAUX_Hist = CoeffDrag * RhoAtm * SrfFlux_Hist.TAUX_Hist
        !!! SrfFlux_Hist.TAUY_Hist = CoeffDrag * RhoAtm * SrfFlux_Hist.TAUY_Hist
    
        END IF
        
    ELSEIF ( Forcing_DT_Unit == 'F' )THEN
    
        !!! CMIP Forcing for spin-up
        IF( (TimTick.Run_AbsIDay <= 365).AND.(TimTick.CMonth == 1).AND.(TimTick.CDay == 1).AND.(mod(TimTick.Run_AbsIT,Day_NT)==1) ) THEN
        
        CountIMon = int4(12);
        StartIMon = int4(12*(TimTick.CYear - Start_CYear) + 1); 
        
        call NC_Read_Mul_IT(Forcing_Path,CMIP_TauX_Forcing_FName,'CMIP_TauX',&
                           SrfFlux_Hist.TAUX_Hist, StartIMon, CountIMon,'N');
        
        call NC_Read_Mul_IT(Forcing_Path,CMIP_TauY_Forcing_FName,'CMIP_TauY',&
                           SrfFlux_Hist.TAUY_Hist, StartIMon, CountIMon,'N');
        
        END IF
        
    ELSEIF ( Forcing_DT_Unit == 'A' )THEN
    
        !!! YYYY.01.01 + first time step, read external forcing data
        IF( (TimTick.CMonth == 1).AND.(TimTick.CDay == 1).AND.(mod(TimTick.Run_AbsIT,Day_NT)==1) ) THEN
    
        CountIDay = int4( 365 + (Leap_Calendar_Onoff) * Math_IsLeapYear(TimTick.CYear) );
        StartIDay = int4( TimTick.Run_AbsIDay ); 
        
        !!! read observational Nino4 & Nino3 daily series
        call NC_Read_Mul_IT(Forcing_Path,OBS_Nudging_FName,'Nino4',&
                           SrfFlux_Hist.Nino4_Hist, StartIDay, CountIDay,'N');
        
        call NC_Read_Mul_IT(Forcing_Path,OBS_Nudging_FName,'Nino3',&
                           SrfFlux_Hist.Nino3_Hist, StartIDay, CountIDay,'N');
        
        
        END IF
    
    ELSEIF ( Forcing_DT_Unit == 'N' )THEN
    
        !!! YYYY.01.01 + first time step, read external forcing data
        IF( (TimTick.CMonth == 1).AND.(TimTick.CDay == 1).AND.(mod(TimTick.Run_AbsIT,Day_NT)==1) ) THEN
    
        CountIDay = int4( 365 + (Leap_Calendar_Onoff) * Math_IsLeapYear(TimTick.CYear) );
        StartIDay = int4( TimTick.Run_AbsIDay ); 
        
        !!! read observational Nino4 & Nino3 daily series
        call NC_Read_Mul_IT(Forcing_Path,OBS_Nudging_FName,'Nino4',&
                           SrfFlux_Hist.Nino4_Hist, StartIDay, CountIDay,'N');
        
        call NC_Read_Mul_IT(Forcing_Path,OBS_Nudging_FName,'Nino3',&
                           SrfFlux_Hist.Nino3_Hist, StartIDay, CountIDay,'N');
        
        
        END IF
        
    END IF !!! Forcing_DT_Unit 
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    IF ( Forcing_DT_Unit == 'C' ) THEN      !!! climtalogical forcing
    
       SrfFlux.TAUX = Clim_OcnStat.TAUXBar (1:GrdInfo.NX+1, 1:GrdInfo.NY,   TimTick.CMonth);
       SrfFlux.TAUY = Clim_OcnStat.TAUYBar (1:GrdInfo.NX,   1:GrdInfo.NY+1, TimTick.CMonth);
      
       !!! SrfFlux.SSH  = SrfFlux_Hist.SSH_Hist(1:GrdInfo.NX+2, 1:GrdInfo.NY+2, TimTick.CMonth);
   
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
    ELSEIF ( Forcing_DT_Unit == 'D' ) THEN  !!! daily forcing
    
       SrfFlux.TAUX = SrfFlux_Hist.TAUX_Hist (1:GrdInfo.NX+1, 1:GrdInfo.NY,   TimTick.YrDaySum);
       SrfFlux.TAUY = SrfFlux_Hist.TAUY_Hist (1:GrdInfo.NX,   1:GrdInfo.NY+1, TimTick.YrDaySum);
       
       IF(Solve_TAUXY_AI_OnOff == 1 ) THEN
          Call Get_TAUXY_AI (GrdInfo, TimTick, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff = 1)
       ELSEIF(Solve_TAUXY_AI_OnOff == 2 ) THEN
          Call Get_TAUXY_AI_EOF (GrdInfo, TimTick, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff = 1)
       ELSE
          Call Get_TAUXY_REG (GrdInfo, TimTick, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff = 1)
       END IF
       
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
       
    ELSEIF ( Forcing_DT_Unit == 'M' ) THEN  !!! monthly forcing
    
       SrfFlux.TAUX = SrfFlux_Hist.TAUX_Hist (1:GrdInfo.NX+1, 1:GrdInfo.NY,   TimTick.CMonth);
       SrfFlux.TAUY = SrfFlux_Hist.TAUY_Hist (1:GrdInfo.NX,   1:GrdInfo.NY+1, TimTick.CMonth);
       
       IF(Solve_TAUXY_AI_OnOff == 1 ) THEN
          Call Get_TAUXY_AI (GrdInfo, TimTick, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff = 1)
       ELSEIF(Solve_TAUXY_AI_OnOff == 2 ) THEN
          Call Get_TAUXY_AI_EOF (GrdInfo, TimTick, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff = 1)
       ELSE
          Call Get_TAUXY_REG (GrdInfo, TimTick, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff = 1) 
       END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    ELSEIF( Forcing_DT_Unit == 'F' ) THEN !!! free coupled run
    
       IF( TimTick.Run_AbsIDay <= 365 ) THEN
       
           !!! clim force spin-up
           !!! SrfFlux.TAUX = Clim_OcnStat.TAUXBar (1:GrdInfo.NX+1, 1:GrdInfo.NY,   TimTick.CMonth);
           !!! SrfFlux.TAUY = Clim_OcnStat.TAUYBar (1:GrdInfo.NX,   1:GrdInfo.NY+1, TimTick.CMonth);
           
           !!! cmip force spin-up
           SrfFlux.TAUX = SrfFlux_Hist.TAUX_Hist (1:GrdInfo.NX+1, 1:GrdInfo.NY,   TimTick.CMonth);
           SrfFlux.TAUY = SrfFlux_Hist.TAUY_Hist (1:GrdInfo.NX,   1:GrdInfo.NY+1, TimTick.CMonth);
           
       
       ELSE !!! TimTick.Run_AbsIDay > 365

          IF(Wind_Nudging_OnOff == 1 ) THEN

            Call Get_TAUXY_REG_NUDGE (GrdInfo, TimTick, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff = 0)

          ELSE !!! Wind_Nudging_OnOff == 0

            Call Get_TAUXY_REG (GrdInfo, TimTick, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff = 0)
          
          END IF
          
       
       END IF !!! TimTick.Run_AbsNday
       
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSEIF( Forcing_DT_Unit == 'A' ) THEN !!! direct assimilation of SSTA
    
      !!! Get Regression SST
      OcnStat.Reg_SST = (Clim_OcnStat.SSTA_Reg_N4) * (SrfFlux_Hist.Nino4_Hist(TimTick.YrDaySum)) + &
                        (Clim_OcnStat.SSTA_Reg_N3) * (SrfFlux_Hist.Nino3_Hist(TimTick.YrDaySum)) + &
                        (Clim_OcnStat.SSTA_Reg_Const)
    
      !!! Fill Galapagos islands (288x64 grid)
      !!! OcnStat.Reg_SST(261:262,32) = (0.5d0) * ( OcnStat.Reg_SST(261:262,31) + OcnStat.Reg_SST(261:262,33) )
    
      OcnStat.Nudging_SST = Nudging_SSTA_Amp * OcnStat.Reg_SST
    
      Call Get_TAUXY_REG (GrdInfo, TimTick, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff = 0)
      
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSEIF( Forcing_DT_Unit == 'N' ) THEN !!! indirect assimilation of SSTA-mapped wind 
    
      OcnStat.Eq_Nudging_T4 = SrfFlux_Hist.Nino4_Hist(TimTick.YrDaySum);
      OcnStat.Eq_Nudging_T3 = SrfFlux_Hist.Nino3_Hist(TimTick.YrDaySum);

      !!! Get Regression SST
      OcnStat.Reg_SST = (Clim_OcnStat.SSTA_Reg_N4) * (SrfFlux_Hist.Nino4_Hist(TimTick.YrDaySum)) + &
                        (Clim_OcnStat.SSTA_Reg_N3) * (SrfFlux_Hist.Nino3_Hist(TimTick.YrDaySum)) + &
                        (Clim_OcnStat.SSTA_Reg_Const)

      !!! Fill Galapagos islands (288x64 grid)
      !!! OcnStat.Reg_SST(261:262,32) = (0.5d0) * ( OcnStat.Reg_SST(261:262,31) + OcnStat.Reg_SST(261:262,33) )

      OcnStat.Nudging_SST = Nudging_SSTA_Amp * OcnStat.Reg_SST

    
      Call Get_TAUXY_REG_NUDGE (GrdInfo, TimTick, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff = 0)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END IF !!! Forcing_DT_Unit 

END SUBROUTINE Ctrl_Add_External_Forcing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
END MODULE MO_PREPOSTPROC
