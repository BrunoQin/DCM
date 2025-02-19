PROGRAM OCEAN_MASTER
    
    use MO_PREPOSTPROC
    use MO_OCEANSOLVER
    use MO_CF23MODEL
    use, intrinsic :: ISO_C_BINDING
    
    IMPLICIT NONE
    
    integer*8 :: SysRun_T0,SysRun_T1
    TYPE(GridInfo) :: GrdInfo
    
    !!!! Tmp_OcnStat for matsuno iteration scheme
    TYPE(OceanStat) :: OcnStat, Tmp_OcnStat 
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(AllFrc_OceanStat) :: AllFrc_OcnStat
    
    TYPE(SurfaceFlux) :: SrfFlux 
    TYPE(SurfaceFlux_Hist) :: SrfFlux_Hist
    
    !!!! Chen-Fang 2023 Simple ENSO ICM param
    TYPE(CFModel_Param) :: CFM_Par
    TYPE(CFModel_Stat)  :: CFM_Stat
    
    integer*8 :: Int_IT, Run_Cycle_ID
    integer*8 :: IM
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! set openmp threads
    call OMP_SET_NUM_THREADS(MY_OMP_THREADS)
    write(*,*),'---------------------------------------------------------------------------'
    write(*,*),'OpenMP THREADS = ',MY_OMP_THREADS
    write(*,*),'---------------------------------------------------------------------------'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
!!! repeatedly run the whole program (initialization + main_cycle)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO Run_Cycle_ID = 1,Run_Cycles 
    
TimTick.Run_Cycle_ID = Run_Cycle_ID
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! initialization
    call Ctrl_Initialization &
    (GrdInfo,Clim_OcnStat,TimTick, OcnStat,Tmp_OcnStat,AllFrc_OcnStat,SrfFlux,SrfFlux_Hist)
    
    call CFModel_Initialization (CFM_Par, CFM_Stat, TimTick)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! assgin inital value
    !DO IM = 1,GrdInfo.NMode
    !    OcnStat.Co_Pres(:,:,IM) = (1.0d0) * Math_TryFunc(GrdInfo.NX+1,GrdInfo.NY+1,'G')
    !END DO 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!! MAIN CYCLE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF(Run_OnOff == 1)then
    
    call system_clock(SysRun_T0)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF(ReRun_OnOff == 1)THEN!!! Read ReStart & ReCover All the Field
       
       call CFModel_Read_ReStart(Restart_Path,CFM_Par,CFM_Stat,TimTick);
    
       !!! ReCover Physical Field (VERY IMPORTANT For the ReStart of the external model ELBOM)
       call CFModel_Recover_Field (TimTick,CFM_Par,CFM_Stat)
       
       IF ( Forcing_DT_Unit == 'F' ) THEN
           call CFModel_Get_Eq_Nudging_Val &
                (GrdInfo,TimTick,CFM_Par,CFM_Stat,Clim_OcnStat,OcnStat) 
       END IF 
       
       call Read_ReStart(ReStart_Path,GrdInfo,OcnStat,TimTick,SrfFlux);
       
       !!! ReCover Surface Field  
       call ReCover_Surface_DynField (GrdInfo,Clim_OcnStat, TimTick, OcnStat, SrfFlux);
       
       !!!  Nino Index
       call Get_NinoIndex (GrdInfo,OcnStat,SrfFlux)
       
    ELSE !!! ReRun_OnOff == 0
     
       call CFModel_ReInit_Stat(CFM_Par,CFM_Stat,TimTick);
       
       !!! ReCover Physical Field (VERY IMPORTANT For the ReStart of the external model ELBOM)
       call CFModel_Recover_Field (TimTick,CFM_Par,CFM_Stat)
       
       IF ( Forcing_DT_Unit == 'F' ) THEN
           call CFModel_Get_Eq_Nudging_Val &
                (GrdInfo,TimTick,CFM_Par,CFM_Stat,Clim_OcnStat,OcnStat) 
       END IF 
        
       call Reinit_Stat(GrdInfo,OcnStat,TimTick,SrfFlux);
       
       !!! ReCover Surface Field  
       call ReCover_Surface_DynField (GrdInfo,Clim_OcnStat, TimTick, OcnStat, SrfFlux);
       
       !!!  Nino Index
       call Get_NinoIndex (GrdInfo,OcnStat,SrfFlux)
       
    END IF !!! ReRun_OnOff
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO Int_IT = (TimTick.Res_AbsIT + 1), (TimTick.Run_AbsNT)
    !!! DO Int_IT = 1,1 
    
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! very important, must start from 1 or (TimTick.Res_AbsIT + 1)
       TimTick.Run_AbsIT = Int_IT; 
       !!! update TimeTick
       call Update_TimeTick (TimTick)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! add external forcing 
       call Ctrl_Add_External_Forcing  &
            (GrdInfo, Clim_OcnStat, OcnStat, TimTick, SrfFlux, SrfFlux_Hist )
       
       !!! add MJO external forcing
       IF ( (Forcing_DT_Unit == 'F').AND.(Add_MJO_Forcing_OnOff == 1) ) THEN
           call Add_MJO_Forcing_TAUXY (GrdInfo,TimTick,OcnStat,Clim_OcnStat,CFM_Stat,SrfFlux)
       END IF
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF (Solver_OnOff == 1 )THEN
       
           !!! The TimeStep Of CF23Model is 12.0d0 hour
           IF ( mod (Int_IT, CFModel_SplitNT ) == 1 ) THEN
              !!! Update CFM_Run_AbsIT
              CFM_Stat.CFM_Run_AbsIT = ceiling (dble(TimTick.Run_AbsIT)/dble(CFModel_SplitNT));
                
              call CFModel_OneStep_Forward (TimTick,CFM_Par,CFM_Stat)
              
              !!! free couple 
              IF ( Forcing_DT_Unit == 'F' ) THEN
                  call CFModel_Get_Eq_Nudging_Val &
                       (GrdInfo,TimTick,CFM_Par,CFM_Stat,Clim_OcnStat,OcnStat) 
              END IF 
              
           END IF !!! mod (Int_IT, CFModel_SplitNT) == 1 

           
           
           IF( CFModel_RunAlone_OnOff == 0 ) THEN
               !!! Ocean Solver
               call OneStep_Forward &
               (GrdInfo,Clim_OcnStat, TimTick, OcnStat, Tmp_OcnStat, AllFrc_OcnStat,SrfFlux)
           END IF !!! CFModel_RunAlone_OnOff
       
       END IF !!! Solver_OnOff
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       !!! Post-processing 
       !!! Monthly Mean 
       IF( MonMe_OutPut_OnOff == 1 ) THEN
       
          call CFModel_Get_MonMe_Stat(TimTick,CFM_Par,CFM_Stat)
       
          IF( CFModel_RunAlone_OnOff == 0 ) then
          call Get_MonMe_Stat(TimTick, GrdInfo, OcnStat, Clim_OcnStat, AllFrc_OcnStat,SrfFlux)
          END IF
          
       END IF
       
       !!! OutPut must before ReStart 
       IF ( OutPut_OnOff == 1) then
       
            call CFModel_Ctrl_OutPut_Stream(TimTick,CFM_Par,CFM_Stat)
            
            IF( CFModel_RunAlone_OnOff == 0 ) then
                call Ctrl_OutPut_Stream (TimTick,GrdInfo,OcnStat,AllFrc_OcnStat, SrfFlux)
            END IF !!! CFModel_RunAlone_OnOff
            
       END IF !!! OutPut_OnOff 
       
       
       IF ( ReStart_OnOff == 1) then
       
            call CFModel_Ctrl_ReStart_Stream (TimTick,CFM_Par,CFM_Stat)
            
            IF( CFModel_RunAlone_OnOff == 0 ) then
                call Ctrl_ReStart_Stream (TimTick,GrdInfo,OcnStat,AllFrc_OcnStat, SrfFlux)
            END IF !!! CFModel_RunAlone_OnOff
            
       END IF !!! ReStart_OnOff 
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    END DO !!! Int_IT
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    call system_clock(SysRun_T1);
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*),'----------------------------------------------------------------------------'
    write(*,*),'CFM_Stat.Eq_T34              = '//trim(Num2Str(CFM_Stat.Eq_T34))
    ! write(*,*),'SrfFlux.TAUX_N4             = '//trim(Num2Str(SrfFlux.TAUX_N4))
    ! write(*,*),'OcnStat.SSTA_N4             = '//trim(Num2Str(OcnStat.SSTA_N4))
    write(*,*),'OcnStat.SSTA_N34             = '//trim(Num2Str(OcnStat.SSTA_N34))
    ! write(*,*),'OcnStat.SSHA_N3             = '//trim(Num2Str(OcnStat.SSHA_N3))
    ! write(*,*),'OcnStat.TSUB_N3             = '//trim(Num2Str(OcnStat.TSUB_N3))
    write(*,*),'----------------------------------------------------------------------------'
    write(*,*),'OpenMP THREADS              = '//trim(Num2Str(MY_OMP_THREADS))
    write(*,*),'TimTick.Run_Cycle_ID        = '//trim(Num2Str(TimTick.Run_Cycle_ID))
    write(*,*),'TimTick.Run_AbsIT           = '//trim(Num2Str(TimTick.Run_AbsIT))
    write(*,*),'MODEL RUNNING SPENDS TIME   = '//trim(Num2Str(dble(SysRun_T1-SysRun_T0)/1.0d6))//' sec'
    write(*,*),'One Run Task Completed'
    write(*,*),'TimTick.Date = '//(TimTick.Date)
    END IF !!! Run_OnOff
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
END DO !!!! Run_Cycle_ID
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
END PROGRAM OCEAN_MASTER
