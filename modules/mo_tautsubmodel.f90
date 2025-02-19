MODULE MO_TAUTSUBMODEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! TAU & TSUB AI-sub model 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE MO_NCTOOLS
USE MO_MATHTOOLS
USE MO_TYPEDEF

USE torch_wrapper

IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION Get_TSUB ( SSHA, SSTA, TCDB) RESULT(TSUB)

    real*8 :: SSHA !!! sea surface height anomaly
    real*8 :: SSTA !!! sea surface temperature anomaly
    real*8 :: TCDB !!! clim thermocline depth
    real*8 :: TSUB
    
    real*8 :: Asub,Bsub
    real*8 :: Hsub,Hstar
    real*8 :: HA,HB,HC,HD
    
     Asub = 11.0d0;
     Bsub= 11.0d0;
      
     Hsub  = 75.0d0;
     Hstar = 60.0d0;
      
     HC = TCDB;
     HD = TCDB;
     HA = (200.0d0) * (SSHA)
      
     TSUB = &
     smheav(+HA) * (+Asub) * ( tanh( (HC+(+HA)-Hsub)/Hstar) - tanh( (HC-Hsub)/Hstar ) ) + &
     smheav(-HA) * (+Bsub) * ( tanh( (HD+(+HA)-Hsub)/Hstar) - tanh( (HD-Hsub)/Hstar ) );

END FUNCTION Get_TSUB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION TSUB_Sigmoid(X) RESULT(Y)
!! thermocline feedback sigmoid function
    
    real*8 :: X,Y
    Y = 1.0d0/(1.0d0 + exp(-0.2d0*X))
    
END FUNCTION TSUB_SigMoid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! PyTorch TSUB Parameterization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE  Get_TSUB_AI(GrdInfo, TimTick, OcnStat, Clim_OcnStat)

    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    
    integer*8,parameter :: AI_PNX = 288, AI_PNY = 64
    
    !!! amplitude modulation factor
    real*8,parameter :: Vanish_LonA = 140.0d0, Vanish_LonB = 180.0d0, Vanish_LonC = 280.0d0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    TYPE(ftorchmodel) :: model;
    INTEGER :: res
    REAL(C_float), dimension(:, :), allocatable :: input
    REAL(C_float) :: output(3, AI_PNY, AI_PNX)
    
    integer*8 :: IX,IY,IEOF
    
    real*8,allocatable :: EOF_TSUB_Coeff(:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( mod (TimTick.Run_AbsIT, (Day_NT)*Tsub_Couple_N_Day) == 1 )THEN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    allocate(EOF_TSUB_Coeff(Tsub_EOF_Trunc_NMode));
    EOF_TSUB_Coeff = 0.0d0;
    
    !!! Assign Value

    allocate (input(AI_PNY, AI_PNX))

    DO IX = 1,AI_PNX
       DO IY = 1,AI_PNY
          input(IY,IX) = REAL ( (TSUB_InPut_Amp) * OcnStat.SSH(IX,IY)  )
       END DO
    END DO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IF ( Write_OutPut_Check == 1 ) THEN
        print *, "Tsub Torch Running"
    END IF
    
    model = TimTick.Tsub_model

    res = tsub_forward(model, input, output)  ! use the model to perform reasoning task

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!! Primitive OutPut
    
    DO IX = 1,AI_PNX
       DO IY = 1,AI_PNY
         
          !!! offline_soda3_tsub.pt:
          !!! (1,IY,IX) : Tmid 
          !!! (2,IY,IX) : Tsrf
          !!! (3,IY,IX) : Tsub
          
          OcnStat.TSUB(IX,IY) = dble( output(2,IY,IX) * GrdInfo.PMask(IX,IY) )

          !!! sign constraint 
          !IF( ( OcnStat.TSUB(IX,IY) * OcnStat.SSH(IX,IY) ) <=0.0d0 ) THEN
          !     OcnStat.TSUB(IX,IY) = (0.0d0) * OcnStat.TSUB(IX,IY)
          !END IF

          !!! OutPut Amp
          OcnStat.TSUB(IX,IY) = (TSUB_OutPut_Amp) * OcnStat.TSUB(IX,IY)
          
       END DO
    END DO
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF( Tsub_EOF_Recon_OnOff == 1) THEN
 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! EOF Reconstruction to remove noise
        !$omp parallel do schedule (static) &
        !$omp default(shared) &
        !$omp private(IEOF)
        DO IEOF = 1,Tsub_EOF_Trunc_NMode
    
           EOF_TSUB_Coeff (IEOF) = sum( OcnStat.TSUB * Clim_OcnStat.EOF_TSUB(1:AI_PNX,1:AI_PNY,IEOF) ); 
    
        END DO
        !$omp end parallel do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        !!! Reconstruction 
        OcnStat.TSUB = 0.0d0;
        DO IEOF = 1,Tsub_EOF_Trunc_NMode
    
           OcnStat.TSUB = OcnStat.TSUB + &
                          EOF_TSUB_Coeff(IEOF) *  Clim_OcnStat.EOF_TSUB(1:AI_PNX,1:AI_PNY,IEOF);
    
        END DO

  
    END IF !!! Tsub_EOF_Recon_Onoff
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF( Tsub_Correct_OnOff == 1) THEN
       
        !!!!!!! Correction Mask (related with mixed ratio based on correlation between TSUB & SSHA)
        DO IX = 1,AI_PNX
            DO IY = 1,AI_PNY

               OcnStat.TSUB(IX,IY) = &
               (        Clim_OcnStat.TSUB_OutPut_Mask(IX,IY,TimTick.CMonth) ) * OcnStat.TSUB(IX,IY) + &
               (1.0d0 - Clim_OcnStat.TSUB_OutPut_Mask(IX,IY,TimTick.CMonth) ) * OcnStat.SST(IX,IY)
 

            END DO
        END DO

    END IF !!! Tsub_Correct_OnOff
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF( Tsub_Slope_Mix_OnOff == 1) THEN
    !!! Mix Ratio with SST
    DO IX = 1,AI_PNX
       DO IY = 1,AI_PNY
        
        !!! Slope Modulation
        IF( ( GrdInfo.PLon(IX)>=Vanish_LonB ).AND. (GrdInfo.PLon(IX)<=Vanish_LonC) ) THEN
          
            OcnStat.TSUB(IX,IY) = &
            ( 1.0d0 - (1.0d0 - Tsub_Slope_Mix_Amp ) * (GrdInfo.PLon(IX)-Vanish_LonC)/(Vanish_LonB-Vanish_LonC) ) * &
            OcnStat.TSUB(IX,IY) + &
            ( (1.0d0 - Tsub_Slope_Mix_Amp ) * (GrdInfo.PLon(IX)-Vanish_LonC)/(Vanish_LonB-Vanish_LonC) ) * &
            OcnStat.SST(IX,IY)
              
          
        ELSEIF( ( GrdInfo.PLon(IX)>= Vanish_LonA).AND. (GrdInfo.PLon(IX)<=Vanish_LonB) ) THEN
          
            OcnStat.TSUB(IX,IY) = &
            ( (1.0d0) -  (1.0d0 - Tsub_Slope_Mix_Amp) * (GrdInfo.PLon(IX) - Vanish_LonA)/(Vanish_LonB-Vanish_LonA) ) * &
            OcnStat.TSUB(IX,IY) + &
            ( (1.0d0 - Tsub_Slope_Mix_Amp) * (GrdInfo.PLon(IX) - Vanish_LonA)/(Vanish_LonB-Vanish_LonA) ) * &
            OcnStat.SST(IX,IY)
          
        ELSE
          
            OcnStat.TSUB(IX,IY) =  OcnStat.TSUB(IX,IY)
          
        ENDIF !!! Slope Modulation
          
       END DO
    END DO
    END IF !!! Tsub_Slope_Mix_OnOff
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! Fill Galapagos islands (288x64 grid)
    OcnStat.TSUB(261:262,32) = (0.5d0) * ( OcnStat.TSUB(261:262,31) + OcnStat.TSUB(261:262,33) )

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END IF !!!! ( mod (TimTick.Run_AbsIT, (Day_NT)*Tsub_N_Day ) == 1 )
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
END SUBROUTINE Get_TSUB_AI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! post process TAUXY 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PostProc_TAUXY (GrdInfo, TimTick, TAUX_TAUY_PGrid, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff)

    TYPE(GridInfo) :: GrdInfo
    TYPE(TimeTick) :: TimTick
    
    real*8,allocatable :: TAUX_TAUY_PGrid(:,:,:) !!! (NX+2,NY+2,2)
    
    TYPE(OceanStat) :: OcnStat
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(SurfaceFlux) :: SrfFlux
    
    integer*8 :: Test_OnOff !!! Offline Test On (1) or Off (0)
    
    integer*8 :: IX,IY,IT
    integer*8 :: NX,NY
    
    !!! sponge mask
    integer*8 :: Tau_Sponge_Num
    integer*8 :: Tau_Sponge_Pow
    real*8,allocatable :: Tau_Sponge_Alp(:) 
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    
    !!!! Set Wind Stress Sponge Layers
    Tau_Sponge_Num = GrdInfo.Sponge_Num
    Tau_Sponge_Pow = GrdInfo.Sponge_Pow
    
    allocate(Tau_Sponge_Alp(Tau_Sponge_Num));
    Tau_Sponge_Alp = GrdInfo.Sponge_Alp
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( Test_OnOff == 1) THEN
    !!!! Test Run
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! PGrid to UGrid, VGrid interpolation
    DO IX = 1,NX+1
        DO IY = 1,NY
        
        SrfFlux.TAUX_Test (IX,IY) = (0.5d0) * (  TAUX_TAUY_PGrid(IX+1,IY+1,1)  + TAUX_TAUY_PGrid(IX,IY+1,1)  );
        
        !!! SrfFlux.TAUX_Test (IX,IY) = &
        !!! ( SrfFlux.TAUX_Test (IX,IY) - Clim_OcnStat.CMIP_TAUXBar (IX,IY,TimTick.CMonth) )
            
          
        END DO
    END DO
        
    DO IX = 1,NX
        DO IY = 1, NY+2
        
        SrfFlux.TAUY_Test (IX,IY) = (0.5d0) * (  TAUX_TAUY_PGrid(IX+1,IY+1,2)  + TAUX_TAUY_PGrid(IX+1,IY,2) );
        
        !!! SrfFlux.TAUY_Test (IX,IY) = &
        !!! ( SrfFlux.TAUY_Test (IX,IY) - Clim_OcnStat.CMIP_TAUYBar (IX,IY,TimTick.CMonth)  )
        
        END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSEIF ( Test_OnOff == 0 ) THEN
    !!! Coupled Run  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO IX = 1,NX+1
        DO IY = 1,NY
            
        SrfFlux.TAUX (IX,IY) = (0.5d0) * ( TAUX_TAUY_PGrid(IX+1,IY+1,1)  + TAUX_TAUY_PGrid(IX,IY+1,1) );
            
        !!!SrfFlux.TAUX (IX,IY) = (Couple_Strength) * &
        !!!( SrfFlux.TAUX (IX,IY) - Clim_OcnStat.CMIP_TAUXBar (IX,IY,TimTick.CMonth) )
            
        SrfFlux.TAUX (IX,IY) = (Couple_Strength) * ( SrfFlux.TAUX(IX,IY) );
            
          
        END DO
    END DO
        
    DO IX = 1,NX
        DO IY = 1,NY+1
        
           
        SrfFlux.TAUY (IX,IY) = (0.5d0) * ( TAUX_TAUY_PGrid(IX+1,IY+1,2)  + TAUX_TAUY_PGrid(IX+1,IY,2) );
           
        !!!SrfFlux.TAUY (IX,IY) = (Couple_Strength) * &
        !!!( SrfFlux.TAUY (IX,IY) - Clim_OcnStat.CMIP_TAUYBar (IX,IY,TimTick.CMonth)  )
           
        SrfFlux.TAUY (IX,IY) = (Couple_Strength) * ( SrfFlux.TAUY(IX,IY) ); 

        
        END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END IF !!! Test_OnOff 
        
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Sponge Mask
    DO IY = 1, Tau_Sponge_Num
        
        !!! TAUX
        SrfFlux.TAUX(1:NX+1,IY) = &
        (1.0d0 - Tau_Sponge_Alp(IY) ) * SrfFlux.TAUX(1:NX+1,IY);
           
        SrfFlux.TAUX(1:NX+1,NY-(IY-1)) = &
        (1.0d0 - Tau_Sponge_Alp(IY) ) * SrfFlux.TAUX(1:NX+1,NY-(IY-1));
           
           
        !!! TAUY 
        SrfFlux.TAUY(1:NX,IY) = &
        (1.0d0 - Tau_Sponge_Alp(IY) ) * SrfFlux.TAUY(1:NX,IY);
           
        SrfFlux.TAUY(1:NX,NY+1-(IY-1)) = &
        (1.0d0 - Tau_Sponge_Alp(IY) ) * SrfFlux.TAUY(1:NX,NY+1-(IY-1));
        
        
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE PostProc_TAUXY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! reconstruct SSH SST by EOF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Recon_SSH_SST ( GrdInfo, TimTick, SSH_SST_PGrid, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff)

    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(SurfaceFlux) :: SrfFlux
    integer*8 :: Test_OnOff !!! Offline Test On (1) or Off (0)
    
    real*8,allocatable :: SSH_SST_PGrid(:,:,:)
    real*8,allocatable :: SSH_SST_PC_Coeff(:,:) 
    real*8,allocatable :: SSH_SST_scale(:)
    
    integer*8 :: IX,IY,IT,IM
    integer*8 :: NX,NY,CMonth,Tau_EOF_NMode
    
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    CMonth = TimTick.CMonth
    Tau_EOF_NMode = Clim_OcnStat.Tau_EOF_NMode
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    allocate(SSH_SST_scale(2));
    SSH_SST_scale = Clim_OcnStat.SSH_SST_scale
    
    allocate(SSH_SST_PC_Coeff(Clim_OcnStat.Tau_EOF_NMode,1) );
    SSH_SST_PC_Coeff = 0.0d0;
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Project onto SSH_SST EOF to get SSH_SST_PC_Coeff
    !!! EOF themselves are unitary vectors due to Climate Data Operator (CDO)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IM)
    DO IM = 1,Tau_EOF_Trunc_NMode
    
       SSH_SST_PC_Coeff(IM,1) = &
       sum( SSH_SST_scale(1) * SSH_SST_PGrid(1:NX+2,1:NY+2,1) * Clim_OcnStat.SSH_SST_EOF(1:NX+2,1:NY+2,1,IM) ) + &
       sum( SSH_SST_scale(2) * SSH_SST_PGrid(1:NX+2,1:NY+2,2) * Clim_OcnStat.SSH_SST_EOF(1:NX+2,1:NY+2,2,IM) );
    
    
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    SSH_SST_PGrid = 0.0d0;
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! reconstrcut SSH SST by EOF
    DO IM = 1, Tau_EOF_Trunc_NMode
       
       !!! SSH
       SSH_SST_PGrid(1:NX+2,1:NY+2,1) = SSH_SST_PGrid(1:NX+2,1:NY+2,1) + &
       SSH_SST_PC_Coeff(IM,1) * (1.0d0/SSH_SST_scale(1)) * Clim_OcnStat.SSH_SST_EOF(1:NX+2,1:NY+2,1,IM)
       
       !!! SST
       SSH_SST_PGrid(1:NX+2,1:NY+2,2) = SSH_SST_PGrid(1:NX+2,1:NY+2,2) + &
       SSH_SST_PC_Coeff(IM,1) * (1.0d0/SSH_SST_scale(2)) * Clim_OcnStat.SSH_SST_EOF(1:NX+2,1:NY+2,2,IM)
    
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE Recon_SSH_SST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! reconstruct TAUX TAUY by EOF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Recon_TAUX_TAUY ( GrdInfo, TimTick, TAUX_TAUY_PGrid, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff)

    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(SurfaceFlux) :: SrfFlux
    integer*8 :: Test_OnOff !!! Offline Test On (1) or Off (0)
    
    real*8,allocatable :: TAUX_TAUY_PGrid(:,:,:)
    real*8,allocatable :: TAUX_TAUY_PC_Coeff(:,:) 
    real*8,allocatable :: TAUX_TAUY_scale(:)
    
    integer*8 :: IX,IY,IT,IM
    integer*8 :: NX,NY,CMonth,Tau_EOF_NMode
    
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    CMonth = TimTick.CMonth
    Tau_EOF_NMode = Clim_OcnStat.Tau_EOF_NMode
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    allocate(TAUX_TAUY_PC_Coeff(Clim_OcnStat.Tau_EOF_NMode,1) );
    TAUX_TAUY_PC_Coeff = 0.0d0;
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Project onto TAUX_TAUY EOF to get TAUX_TAUY_PC_Coeff
    !!! EOF themselves are unitary vectors due to Climate Data Operator (CDO)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IM)
    DO IM = 1,Tau_EOF_Trunc_NMode
    
       TAUX_TAUY_PC_Coeff(IM,1) = &
       sum( TAUX_TAUY_PGrid(1:NX+2,1:NY+2,1) * Clim_OcnStat.TAUX_TAUY_EOF(1:NX+2,1:NY+2,1,IM) ) + &
       sum( TAUX_TAUY_PGrid(1:NX+2,1:NY+2,2) * Clim_OcnStat.TAUX_TAUY_EOF(1:NX+2,1:NY+2,2,IM) );
    
    
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    TAUX_TAUY_PGrid = 0.0d0;
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! reconstrcut TAUX TAUY by EOF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO IM = 1, Tau_EOF_Trunc_NMode
       
       !!! TAUX
       TAUX_TAUY_PGrid(1:NX+2,1:NY+2,1) = TAUX_TAUY_PGrid(1:NX+2,1:NY+2,1) + &
       TAUX_TAUY_PC_Coeff(IM,1) * Clim_OcnStat.TAUX_TAUY_EOF(1:NX+2,1:NY+2,1,IM)
       
       !!! TAUY
       TAUX_TAUY_PGrid(1:NX+2,1:NY+2,2) = TAUX_TAUY_PGrid(1:NX+2,1:NY+2,2) + &
       TAUX_TAUY_PC_Coeff(IM,1) * Clim_OcnStat.TAUX_TAUY_EOF(1:NX+2,1:NY+2,2,IM)
    
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE Recon_TAUX_TAUY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Traditional EOF Linear regression model to get TAUXY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Get_TAUXY_EOF (GrdInfo, TimTick, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff)

    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(SurfaceFlux) :: SrfFlux
    integer*8 :: Test_OnOff !!! Offline Test On (1) or Off (0)
    
    integer*8 :: IX,IY,IT,IM
    integer*8 :: NX,NY,CMonth,Tau_EOF_NMode
    
    real*8,allocatable :: TAUX_TAUY_PGrid(:,:,:)
    
    !!! they are 2-D matrix (Tau_EOF_NMode,1) for matmul operation
    real*8,allocatable :: SSH_SST_PC_Coeff(:,:) 
    real*8,allocatable :: TAUX_TAUY_PC_Coeff(:,:)
    
    real*8,allocatable :: SSH_SST_scale(:)
    
    real*8,allocatable :: Tau_AMat(:,:),Tau_BVec(:,:)
    
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    CMonth = TimTick.CMonth
    Tau_EOF_NMode = Clim_OcnStat.Tau_EOF_NMode

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF ( mod (TimTick.Run_AbsIT, (Day_NT)*Tau_Couple_N_Day) == 1 )THEN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! allocate temporary computation space
    
    allocate(TAUX_TAUY_PGrid(NX+2,NY+2,2));
    TAUX_TAUY_PGrid = 0.0d0;
    
    allocate(SSH_SST_PC_Coeff  (Clim_OcnStat.Tau_EOF_NMode,1) );
    allocate(TAUX_TAUY_PC_Coeff(Clim_OcnStat.Tau_EOF_NMode,1) );
    
    SSH_SST_PC_Coeff = 0.0d0;
    TAUX_TAUY_PC_Coeff = 0.0d0;
    
    allocate(SSH_SST_scale(2));
    SSH_SST_scale = Clim_OcnStat.SSH_SST_scale
    
    allocate(Tau_AMat(Tau_EOF_NMode,Tau_EOF_NMode));
    allocate(Tau_BVec(Tau_EOF_NMode,1));
    
    Tau_AMat(1:Tau_EOF_NMode,1:Tau_EOF_NMode) = &
    Clim_OcnStat.Tau_AMat_12Mon(1:Tau_EOF_NMode,1:Tau_EOF_NMode,CMonth)
    
    Tau_BVec(1:Tau_EOF_NMode,1) = Clim_OcnStat.Tau_BVec_12Mon(1:Tau_EOF_NMode,CMonth)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Project onto SSH_SST EOF to get SSH_SST_PC_Coeff
    !!! EOF themselves are unitary vectors due to Climate Data Operator (CDO)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IM)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO IM = 1,Tau_EOF_Trunc_NMode
    
       SSH_SST_PC_Coeff(IM,1) = &
       sum( SSH_SST_scale(1) * OcnStat.SSH(1:NX+2,1:NY+2) * Clim_OcnStat.SSH_SST_EOF(1:NX+2,1:NY+2,1,IM) ) + &
       sum( SSH_SST_scale(2) * OcnStat.SST(1:NX+2,1:NY+2) * Clim_OcnStat.SSH_SST_EOF(1:NX+2,1:NY+2,2,IM) );
    
    
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Use linear regression model to get TAUX_TAUY_PC_Coeff
    TAUX_TAUY_PC_Coeff = matmul(Tau_AMat, SSH_SST_PC_Coeff ) + Tau_BVec
    print*,'Linear EOF Tau Model Running'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! reconstruction according to TAUX_TAUY_PC_Coeff
    DO IM = 1, Tau_EOF_Trunc_NMode
       
       !!! TAUX
       TAUX_TAUY_PGrid(1:NX+2,1:NY+2,1) = TAUX_TAUY_PGrid(1:NX+2,1:NY+2,1) + &
       TAUX_TAUY_PC_Coeff(IM,1) * Clim_OcnStat.TAUX_TAUY_EOF(1:NX+2,1:NY+2,1,IM)
       
       !!! TAUY
       TAUX_TAUY_PGrid(1:NX+2,1:NY+2,2) = TAUX_TAUY_PGrid(1:NX+2,1:NY+2,2) + &
       TAUX_TAUY_PC_Coeff(IM,1) * Clim_OcnStat.TAUX_TAUY_EOF(1:NX+2,1:NY+2,2,IM)
    
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! post proc TAUX_TAUY_PGrid to SrfFlux.TAUX & SrfFlux.TAUY
    CALL PostProc_TAUXY (GrdInfo, TimTick, TAUX_TAUY_PGrid, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END IF !!!! IF ( mod (TimTick.Run_AbsIT, (Day_NT)*Tau_Couple_N_Day) == 1 )THEN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END SUBROUTINE Get_TAUXY_EOF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  PC2PC EOF AI model to get TAUXY 
!!!  MAX Tau_EOF_Trunc_NMode = 64
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Get_TAUXY_AI_EOF (GrdInfo, TimTick, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff)

    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(SurfaceFlux) :: SrfFlux
    integer*8 :: Test_OnOff !!! Offline Test On (1) or Off (0)
    
    integer*8 :: IX,IY,IT,IM
    integer*8 :: NX,NY,CMonth,Tau_EOF_NMode
    
    real*8,allocatable :: TAUX_TAUY_PGrid(:,:,:)
    
    !!! they are 2-D matrix (Tau_EOF_NMode,1) for matmul operation
    real*8,allocatable :: SSH_SST_PC_Coeff(:,:) 
    real*8,allocatable :: TAUX_TAUY_PC_Coeff(:,:)
    real*8,allocatable :: SSH_SST_scale(:)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer*8,parameter :: AI_PNX = 288, AI_PNY = 64
    REAL(C_float), dimension(:, :, :), allocatable :: input
    REAL(C_float) :: output(2, AI_PNY, AI_PNX)
    TYPE(ftorchmodel) :: model;
    INTEGER :: res
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    CMonth = TimTick.CMonth
    Tau_EOF_NMode = Clim_OcnStat.Tau_EOF_NMode

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF ( mod (TimTick.Run_AbsIT, (Day_NT)*Tau_Couple_N_Day) == 1 )THEN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! allocate temporary computation space
    allocate( input(2, AI_PNY, AI_PNX ) );
    input = REAL(0.0d0);
    
    allocate(TAUX_TAUY_PGrid(NX+2,NY+2,2));
    TAUX_TAUY_PGrid = 0.0d0;
    
    allocate(SSH_SST_PC_Coeff  (Clim_OcnStat.Tau_EOF_NMode,1) );
    allocate(TAUX_TAUY_PC_Coeff(Clim_OcnStat.Tau_EOF_NMode,1) );
    SSH_SST_PC_Coeff = 0.0d0;
    TAUX_TAUY_PC_Coeff = 0.0d0;
    
    allocate(SSH_SST_scale(2));
    SSH_SST_scale = Clim_OcnStat.SSH_SST_scale
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Project onto SSH_SST EOF to get SSH_SST_PC_Coeff
    !!! EOF themselves are unitary vectors due to Climate Data Operator (CDO)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IM)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO IM = 1,Tau_EOF_Trunc_NMode
    
       SSH_SST_PC_Coeff(IM,1) = &
       sum( SSH_SST_scale(1) * OcnStat.SSH(1:NX+2,1:NY+2) * Clim_OcnStat.SSH_SST_EOF(1:NX+2,1:NY+2,1,IM) ) + &
       sum( SSH_SST_scale(2) * OcnStat.SST(1:NX+2,1:NY+2) * Clim_OcnStat.SSH_SST_EOF(1:NX+2,1:NY+2,2,IM) );
    
    
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    input(1,1:AI_PNY,1) = REAL(SSH_SST_PC_Coeff(1:AI_PNY,1))
    
    print *, "EOF Tau Torch Running"
    model = TimTick.Tau_model
    res = tau_forward(model, input, output)  ! use the model to perform reasoning task
    
    TAUX_TAUY_PC_Coeff(1:AI_PNY,1) = dble( output(1,1:AI_PNY,1) );
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! reconstruction according to TAUX_TAUY_PC_Coeff
    DO IM = 1, Tau_EOF_Trunc_NMode
       
       !!! TAUX
       TAUX_TAUY_PGrid(1:NX+2,1:NY+2,1) = TAUX_TAUY_PGrid(1:NX+2,1:NY+2,1) + &
       TAUX_TAUY_PC_Coeff(IM,1) * Clim_OcnStat.TAUX_TAUY_EOF(1:NX+2,1:NY+2,1,IM)
       
       !!! TAUY
       TAUX_TAUY_PGrid(1:NX+2,1:NY+2,2) = TAUX_TAUY_PGrid(1:NX+2,1:NY+2,2) + &
       TAUX_TAUY_PC_Coeff(IM,1) * Clim_OcnStat.TAUX_TAUY_EOF(1:NX+2,1:NY+2,2,IM)
    
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! post proc TAUX_TAUY_PGrid to SrfFlux.TAUX & SrfFlux.TAUY
    CALL PostProc_TAUXY (GrdInfo, TimTick, TAUX_TAUY_PGrid, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END IF !!!! IF ( mod (TimTick.Run_AbsIT, (Day_NT)*Tau_Couple_N_Day) == 1 )THEN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END SUBROUTINE Get_TAUXY_AI_EOF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Get_TAUXY_AI( GrdInfo, TimTick, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff)

    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(SurfaceFlux) :: SrfFlux
    integer*8 :: Test_OnOff !!! Offline Test On (1) or Off (0)
    
    integer*8,parameter :: AI_PNX = 288, AI_PNY = 64
    integer*8,parameter :: NX = 286, NY = 62
    
    REAL(C_float), dimension(:, :, :), allocatable :: input
    REAL(C_float) :: output(2, AI_PNY, AI_PNX)
    
    real*8,allocatable :: SSH_SST_PGrid(:,:,:)
    real*8,allocatable :: TAUX_TAUY_PGrid(:,:,:)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    TYPE(ftorchmodel) :: model;
    INTEGER :: res
    integer*8 :: IX,IY,IT

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( mod (TimTick.Run_AbsIT, (Day_NT)*Tau_Couple_N_Day) == 1 )THEN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        !!! Reconstruct SSH SST by EOF
        allocate(SSH_SST_PGrid(AI_PNX,AI_PNY,2));

        SSH_SST_PGrid(          1 : AI_PNX/2, 1:AI_PNY, 1) = OcnStat.SSTA_N4
        
        SSH_SST_PGrid( AI_PNX/2+1 : AI_PNX  , 1:AI_PNY, 1) = OcnStat.SSTA_N3

        SSH_SST_PGrid(          1 : AI_PNX  , 1:AI_PNY, 2) = dble(TimTick.CMonth)
        
        !!! CALL Recon_SSH_SST ( GrdInfo, TimTick, SSH_SST_PGrid, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff)
        
        !!! write(*,*),'OcnStat.SSTA_N4 = ',OcnStat.SSTA_N4
        !!! write(*,*),'OcnStat.SSTA_N3 = ',OcnStat.SSTA_N3

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate( input(2, AI_PNY, AI_PNX ) );
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !$omp parallel do schedule (static) &
        !$omp default(shared) &
        !$omp private(IX,IY)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO IX = 1,AI_PNX
           DO IY = 1,AI_PNY
           
              !!! input(1, IY, IX ) = REAL( OcnStat.SSH ( IX, IY ) + Clim_OcnStat.CMIP_SSHBar (IX,IY,TimTick.CMonth) );
              !!! input(2, IY, IX ) = REAL( OcnStat.SST ( IX, IY ) + Clim_OcnStat.CMIP_SSTBar (IX,IY,TimTick.CMonth) );
              !!! input(3, IY, IX ) = REAL( TimTick.CMonth );
              
              input(1, IY, IX ) = REAL( SSH_SST_PGrid ( IX, IY, 1)  );
              input(2, IY, IX ) = REAL( SSH_SST_PGrid ( IX, IY, 2)  );
    
           END DO
        END DO
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !$omp end parallel do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        print *, "Tau Torch Running"
        model = TimTick.Tau_model

        res = tau_forward(model, input, output)  ! use the model to perform reasoning task
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
        !!! allocate double tauxy output
        allocate(TAUX_TAUY_PGrid(NX+2,NY+2,2));
        TAUX_TAUY_PGrid = 0.0d0;
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !$omp parallel do schedule (static) &
        !$omp default(shared) &
        !$omp private(IX,IY)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! Tau_PMask + Tau_Bias
        DO IX = 1,NX+2
            DO IY = 1,NY+2
           
            !!! successful operation
            !!! TAUX_TAUY_PGrid(IX,IY,1) =  dble( output(1,IY,IX) * REAL(GrdInfo.PMask(IX,IY)) + &
            !!!                              REAL(Couple_TAUX_Bias) * REAL(GrdInfo.PMask(IX,IY)) );
            !!!
            !!! TAUX_TAUY_PGrid(IX,IY,2) =  dble( output(2,IY,IX) * REAL(GrdInfo.PMask(IX,IY)) );
            
            TAUX_TAUY_PGrid(IX,IY,1) =  dble( output(1,IY,IX) );
            
            TAUX_TAUY_PGrid(IX,IY,2) =  dble( output(2,IY,IX) );
           
            END DO
        END DO
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !$omp end parallel do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        !!! Reconstruct TAUXY by EOF
        !!! CALL Recon_TAUX_TAUY ( GrdInfo, TimTick, TAUX_TAUY_PGrid, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff)
        
        
        !!! post proc TAUX_TAUY_PGrid to SrfFlux.TAUX & SrfFlux.TAUY
        CALL PostProc_TAUXY (GrdInfo, TimTick, TAUX_TAUY_PGrid, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff)
        
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END IF !!! IF ( mod (TimTick.Run_AbsIT, (Day_NT)*Tau_Couple_N_Day) == 1 )THEN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

END SUBROUTINE Get_TAUXY_AI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Get_TAUXY_REG( GrdInfo, TimTick, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff)

!!! Regression Model onto Nino4 & Nino3 Index

    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(SurfaceFlux) :: SrfFlux
    integer*8 :: Test_OnOff !!! Offline Test On (1) or Off (0)
    
    integer*8,parameter :: AI_PNX = 288, AI_PNY = 64
    !!! integer*8,parameter :: NX = 286, NY = 62
    
    
    integer*8 :: IX,IY,IT
    real*8,allocatable :: TAUX_TAUY_PGrid(:,:,:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( mod (TimTick.Run_AbsIT, (Day_NT)*Tau_Couple_N_Day) == 1 )THEN
    
    
        !!! allocate double tauxy output
        allocate(TAUX_TAUY_PGrid(AI_PNX,AI_PNY,2));
        TAUX_TAUY_PGrid = 0.0d0;
    
    
        TAUX_TAUY_PGrid(1:AI_PNX,1:AI_PNY,1) = &
        Clim_OcnStat.TAUX_Reg_N4(1:AI_PNX,1:AI_PNY,TimTick.CMonth) * (OcnStat.SSTA_N4 + Reg_N4_Bias) + &
        Clim_OcnStat.TAUX_Reg_N3(1:AI_PNX,1:AI_PNY,TimTick.CMonth) * (OcnStat.SSTA_N3 + Reg_N3_Bias) + &
        Clim_OcnStat.TAUX_Reg_Const(1:AI_PNX,1:AI_PNY,TimTick.CMonth) 
        
        TAUX_TAUY_PGrid(1:AI_PNX,1:AI_PNY,1) = TAUX_TAUY_PGrid(1:AI_PNX,1:AI_PNY,1) + Reg_Couple_TAUX_Bias;
        
    
        TAUX_TAUY_PGrid(1:AI_PNX,1:AI_PNY,2) = &
        Clim_OcnStat.TAUY_Reg_N4(1:AI_PNX,1:AI_PNY,TimTick.CMonth) * (OcnStat.SSTA_N4 + Reg_N4_Bias) + &
        Clim_OcnStat.TAUY_Reg_N3(1:AI_PNX,1:AI_PNY,TimTick.CMonth) * (OcnStat.SSTA_N3 + Reg_N3_Bias) + &
        Clim_OcnStat.TAUY_Reg_Const(1:AI_PNX,1:AI_PNY,TimTick.CMonth)
        
    
        !!! post proc TAUX_TAUY_PGrid to SrfFlux.TAUX & SrfFlux.TAUY
        CALL PostProc_TAUXY (GrdInfo, TimTick, TAUX_TAUY_PGrid, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff)
    
        !!! deallocate TAUX_TAUY_PGrid
        deallocate(TAUX_TAUY_PGrid)
    
    
    END IF !!! IF ( mod (TimTick.Run_AbsIT, (Day_NT)*Tau_Couple_N_Day) == 1 )THEN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
END SUBROUTINE Get_TAUXY_REG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Get_TAUXY_REG_NUDGE( GrdInfo, TimTick, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff)

!!! Regression Model onto Nino4 & Nino3 Index

    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(SurfaceFlux) :: SrfFlux
    integer*8 :: Test_OnOff !!! Offline Test On (1) or Off (0)
    
    integer*8,parameter :: AI_PNX = 288, AI_PNY = 64
    !!! integer*8,parameter :: NX = 286, NY = 62
    
    
    integer*8 :: IX,IY,IT
    
    real*8,allocatable :: Ext_TAUX_TAUY_PGrid(:,:,:) !!! External Forcing From CF23
    real*8,allocatable :: Cpl_TAUX_TAUY_PGrid(:,:,:) !!! Coupled Forcing From ELBOM (ENBOM)
    real*8,allocatable :: TAUX_TAUY_PGrid(:,:,:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( mod (TimTick.Run_AbsIT, (Day_NT)*Tau_Couple_N_Day) == 1 )THEN
    
        !!! allocate double tauxy output
        allocate(Ext_TAUX_TAUY_PGrid(AI_PNX, AI_PNY, 2));
        allocate(Cpl_TAUX_TAUY_PGrid(AI_PNX, AI_PNY, 2));
        allocate(   TAUX_TAUY_PGrid (AI_PNX, AI_PNY, 2));      
        
        Ext_TAUX_TAUY_PGrid = 0.0d0;
        Cpl_TAUX_TAUY_PGrid = 0.0d0;
        TAUX_TAUY_PGrid = 0.0d0;
        
        call  Compute_TAUX_TAUY_PGrid &
              (Ext_TAUX_TAUY_PGrid, Clim_OcnStat, OcnStat.Eq_Nudging_T4, OcnStat.Eq_Nudging_T3, TimTick.CMonth)
        
        call  Compute_TAUX_TAUY_PGrid &
              (Cpl_TAUX_TAUY_PGrid, Clim_OcnStat, OcnStat.SSTA_N4, OcnStat.SSTA_N3, TimTick.CMonth)
        
        !!! very important, mixing External Wind + Internal Wind
        TAUX_TAUY_PGrid = (Wind_Nudging_Ratio) * Ext_TAUX_TAUY_PGrid + &
                          (1.0d0 - Wind_Nudging_Ratio) * Cpl_TAUX_TAUY_PGrid
    
        !!! post proc TAUX_TAUY_PGrid to SrfFlux.TAUX & SrfFlux.TAUY
        CALL PostProc_TAUXY (GrdInfo, TimTick, TAUX_TAUY_PGrid, OcnStat, Clim_OcnStat, SrfFlux, Test_OnOff)
    
        !!! deallocate TAUX_TAUY_PGrid
        deallocate(Ext_TAUX_TAUY_PGrid, Cpl_TAUX_TAUY_PGrid, TAUX_TAUY_PGrid)
    
    
    END IF !!! IF ( mod (TimTick.Run_AbsIT, (Day_NT)*Tau_Couple_N_Day) == 1 )THEN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
END SUBROUTINE Get_TAUXY_REG_NUDGE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Compute_TAUX_TAUY_PGrid &
           (TAUX_TAUY_PGrid, Clim_OcnStat, SSTA_N4, SSTA_N3, CMonth)

    real*8,allocatable :: TAUX_TAUY_PGrid(:,:,:)
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    
    integer*8 :: AI_PNX, AI_PNY, CMonth
    real*8 :: SSTA_N4, SSTA_N3
    
    AI_PNX = size(TAUX_TAUY_PGrid,1);
    AI_PNY = size(TAUX_TAUY_PGrid,2);
    
    !!! TAUX_TAUY_PGrid = 0.0d0; !!! for safety
    
    
    TAUX_TAUY_PGrid(1:AI_PNX,1:AI_PNY,1) = &
    Clim_OcnStat.TAUX_Reg_N4(1:AI_PNX,1:AI_PNY,CMonth) * (SSTA_N4 + Reg_N4_Bias) + &
    Clim_OcnStat.TAUX_Reg_N3(1:AI_PNX,1:AI_PNY,CMonth) * (SSTA_N3 + Reg_N3_Bias) + &
    Clim_OcnStat.TAUX_Reg_Const(1:AI_PNX,1:AI_PNY,CMonth) 
        
    TAUX_TAUY_PGrid(1:AI_PNX,1:AI_PNY,1) = TAUX_TAUY_PGrid(1:AI_PNX,1:AI_PNY,1) + Reg_Couple_TAUX_Bias;
        
    TAUX_TAUY_PGrid(1:AI_PNX,1:AI_PNY,2) = &
    Clim_OcnStat.TAUY_Reg_N4(1:AI_PNX,1:AI_PNY,CMonth) * (SSTA_N4 + Reg_N4_Bias) + &
    Clim_OcnStat.TAUY_Reg_N3(1:AI_PNX,1:AI_PNY,CMonth) * (SSTA_N3 + Reg_N3_Bias) + &
    Clim_OcnStat.TAUY_Reg_Const(1:AI_PNX,1:AI_PNY,CMonth)
    

END SUBROUTINE Compute_TAUX_TAUY_PGrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Recon_MJO_Forcing_PsTauXY( GrdInfo, TimTick, OcnStat, Clim_OcnStat, SrfFlux, YearDay)

    !!!! reconstruct MJO external forcing

    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(SurfaceFlux) :: SrfFlux
    
    integer*8 :: YearDay !!! from 1 to 365
    
    integer*8 :: IX,IY,IZ,IT,IEEOF
    integer*8 :: NX,NY,NEEOF,CFM_NOcn

    !!! west east phase assymetry
    !!! real*8,parameter :: MJO_Asm_Epi = 1.0d-8 
    !!! real*8,parameter :: MJO_Asm_Fac = 0.95d0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    
    CFM_NOcn = SrfFlux.CFM_NOcn
    NEEOF    = SrfFlux.NEEOF
    
    !!! zero for safety
    SrfFlux.MJO_PsTauX = 0.0d0;
    SrfFlux.MJO_PsTauY = 0.0d0;
    SrfFlux.MJO_KelCoeff   = 0.0d0;

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel sections
    !$omp section
    DO IEEOF = 1,SrfFlux.NEEOF
    
        SrfFlux.MJO_PsTauX = SrfFlux.MJO_PsTauX + ( SrfFlux.PC_MJO_PsTauXY(IEEOF) ) * (0.5d0) * &
        ( Clim_OcnStat.EEOF_PsTauX(1:NX+1,2:NY+1,YearDay,IEEOF) + Clim_OcnStat.EEOF_PsTauX(2:NX+2,2:NY+1,YearDay,IEEOF) );
        
    END DO

    !$omp section
    DO IEEOF = 1,SrfFlux.NEEOF
        
        SrfFlux.MJO_PsTauY = SrfFlux.MJO_PsTauY + ( SrfFlux.PC_MJO_PsTauXY(IEEOF) ) * (0.5d0) * &
        ( Clim_OcnStat.EEOF_PsTauY(2:NX+1,1:NY+1,YearDay,IEEOF) + Clim_OcnStat.EEOF_PsTauY(2:NX+1,2:NY+2,YearDay,IEEOF) );
        
    END DO
    
    !$omp end parallel sections
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    DO IEEOF = 1,SrfFlux.NEEOF
    
    
       SrfFlux.MJO_KelCoeff = SrfFlux.MJO_KelCoeff + SrfFlux.PC_MJO_PsTauXY(IEEOF) * &
                              ( Clim_OcnStat.EEOF_KelCoeff(1:CFM_NOcn,YearDay,IEEOF) )
    
    
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! gain positive value
    !!!
    !!! SrfFlux.MJO_KelCoeff = &
    !!! (            0.5d0) * ( SrfFlux.MJO_KelCoeff + sqrt(SrfFlux.MJO_KelCoeff**(2.0d0) + MJO_Asm_Epi) ) + &
    !!! (MJO_Asm_Fac*0.5d0) * ( SrfFlux.MJO_KelCoeff - sqrt(SrfFlux.MJO_KelCoeff**(2.0d0) + MJO_Asm_Epi) );

    SrfFlux.MJO_KelCoeff = (1.0d0)*(0.5d0)*(SrfFlux.MJO_KelCoeff + abs(SrfFlux.MJO_KelCoeff)) + &
                           (Add_MJO_EastWind_Ratio)*(0.5d0)*(SrfFlux.MJO_KelCoeff - abs(SrfFlux.MJO_KelCoeff))
   
    !!! add bias
    SrfFlux.MJO_KelCoeff = SrfFlux.MJO_KelCoeff + (MJO_WestWind_Bias/(CoeffDrag * RhoAtm))*(PsTauXY_To_KelCoeff)
    
    SrfFlux.MJO_PsTauX   = SrfFlux.MJO_PsTauX + (MJO_WestWind_Bias/(CoeffDrag*RhoAtm))


END SUBROUTINE Recon_MJO_Forcing_PsTauXY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Add_MJO_Forcing_TAUXY ( GrdInfo, TimTick, OcnStat, Clim_OcnStat, CFM_Stat, SrfFlux)

    !!!! reconstruct MJO external forcing

    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(SurfaceFlux) :: SrfFlux
    TYPE(CFModel_Stat)  :: CFM_Stat
    
    integer*8 :: YearDay !!! from 1 to 365
    
    integer*8 :: IX,IY,IZ,IT,IEEOF
    integer*8 :: NX,NY,NEEOF,CFM_NOcn
    
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    
    CFM_NOcn = SrfFlux.CFM_NOcn
    NEEOF    = SrfFlux.NEEOF
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (TimTick.Run_Cycle_ID == 1) THEN
        SrfFlux.PC_MJO_PsTauXY = 0.0d0;
    
    ELSEIF( TimTick.Run_Cycle_ID >= 2) THEN
    
        SrfFlux.PC_MJO_PsTauXY = Clim_OcnStat.PC_PsTauXY(1:NEEOF,TimTick.Run_Cycle_ID-1);
        
        IF( (TimTick.Run_AbsIT - TimTick.Res_AbsIT) > (Day_NT)*180 ) THEN
             SrfFlux.PC_MJO_PsTauXY = 0.0d0;
        END IF
        
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! Must Add Forcing From 01.01 to 6.30 or 12.31
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF ( TimTick.Run_Cycle_ID >= 2) THEN
    IF ( mod (TimTick.Run_AbsIT, Day_NT) == 1 )THEN
    
         YearDay = mod( TimTick.Run_AbsIDay-1, 365) + 1;
         
         !!!! reconstruct MJO-type wind stress anomalies 
         call Recon_MJO_Forcing_PsTauXY( GrdInfo, TimTick, OcnStat, Clim_OcnStat, SrfFlux, YearDay)

         IF( (TimTick.Run_AbsIT - TimTick.Res_AbsIT) > (Day_NT)*365 ) THEN 
             SrfFlux.MJO_PsTauX   = 0.0d0;
             SrfFlux.MJO_PsTauY   = 0.0d0;
             SrfFlux.MJO_KelCoeff = 0.0d0;
         END IF
         
         !!! Add Forcing

         IF(Add_ENBOM_MJO_Forcing_OnOff == 1 )THEN
            SrfFlux.TAUX = SrfFlux.TAUX + CoeffDrag * RhoAtm * SrfFlux.MJO_PsTauX 
            !!! SrfFlux.TAUY = SrfFlux.TAUY + CoeffDrag * RhoAtm * SrfFlux.MJO_PsTauY
         ENDIF


         IF(Add_CF23_MJO_Forcing_OnOff  == 1 ) THEN
            CFM_Stat.Ex_Ocn_Str = CoeffDrag * RhoAtm * SrfFlux.MJO_KelCoeff
         ENDIF

    
    END IF !!! mod (TimTick.Run_AbsIT, Day_NT) == 1
    END IF !!! TimTick.Run_Cycle_ID >= 2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

END SUBROUTINE Add_MJO_Forcing_TAUXY 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE MO_TAUTSUBMODEL
