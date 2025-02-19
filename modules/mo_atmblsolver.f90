MODULE MO_ATMBLSOLVER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Gill-type atmospheric boundary layer solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE MO_NCTOOLS
USE MO_MATHTOOLS
USE MO_TYPEDEF
USE MO_TAUTSUBMODEL

IMPLICIT NONE

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Get_AllFrc_SurfaceFlux  &
           (GrdInfo,Clim_OcnStat, TimTick, SrfFlux, AllFrc_SrfFlux, OcnStat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! ArakaWa C Grid, Equi-Distant Nodes 
!!!
!!!                                        U(IX-1,IY+0)
!!!                           V(IX-1,IY+0) P(IX+0,IY+1) V(IX-1,IY+1)
!!!              U(IX+0,IY-1)              U(IX+0,IY+0)              U(IX+0,IY+1)
!!! V(IX+0,IY-1) P(IX+1,IY+0) V(IX+0,IY+0) P(IX+1,IY+1) V(IX+0,IY+1) P(IX+1,IY+2) V(IX+0,IY+2)
!!!              U(IX+1,IY-1)              U(IX+1,IY+0)              U(IX+1,IY+1)
!!!                           V(IX+1,IY+0) P(IX+2,IY+1) V(IX+1,IY+1) 
!!!                                        U(IX+2,IY+0)                      
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(AllFrc_SurfaceFlux) :: AllFrc_SrfFlux
    TYPE(SurfaceFlux) :: SrfFlux
    
    integer*8 :: NMode
    integer*8 :: IM !!! indicate the index of baroclinic mode
    
    integer*8 :: NX,NY,NZ
    integer*8 :: IX,IY,IZ
    
    integer*8 :: CMonth !!! the calendar month
    
    real*8 :: DX,DY
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    NZ = GrdInfo.NZ
    NMode = GrdInfo.NMode
    
    DX = GrdInfo.ReDX
    DY = GrdInfo.ReDY
    
    CMonth = TimTick.CMonth
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!$omp parallel do schedule (static) &
    !!$omp default(shared) &
    !!$omp private(IX,IY)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! F(Uwnd)
    DO IX = 2,NX
       DO IY = 1,NY
       
        !!! pressure gradient force
        AllFrc_SrfFlux.AllFrc_Gill_Uwnd(IX,IY) = &
        (-1.d0/DX) * (SrfFlux.Gill_Pres(IX+1,IY+1) - SrfFlux.Gill_Pres(IX,IY+1) );
          
        !!!! Coriolis force
        AllFrc_SrfFlux.AllFrc_Gill_Uwnd(IX,IY) = AllFrc_SrfFlux.AllFrc_Gill_Uwnd(IX,IY) + &
        ( EarthBeta * GrdInfo.UReY(IY) ) * (0.25d0) * &
        ( SrfFlux.Gill_Vwnd(IX-1,IY+0) + SrfFlux.Gill_Vwnd(IX-1,IY+1) + &
          SrfFlux.Gill_Vwnd(IX+0,IY+0) + SrfFlux.Gill_Vwnd(IX+0,IY+1) );
        
        !!! Uwnd Damping
        AllFrc_SrfFlux.AllFrc_Gill_Uwnd(IX,IY) = AllFrc_SrfFlux.AllFrc_Gill_Uwnd(IX,IY) + &
        (-1.0d0) * (Gill_UVwnd_Damp) * SrfFlux.Gill_Uwnd(IX,IY)
       
       END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!$omp parallel do schedule (static) &
    !!$omp default(shared) &
    !!$omp private(IX,IY)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! F(Vwnd)
    DO IX = 1,NX
       DO IY = 2,NY
       
          !!! pressure gradient force
          AllFrc_SrfFlux.AllFrc_Gill_Vwnd(IX,IY) = &
          (-1.d0/DY) * (SrfFlux.Gill_Pres(IX+1,IY+1) - SrfFlux.Gill_Pres(IX+1,IY) );
          
          !!!!! Coriolis force
          AllFrc_SrfFlux.AllFrc_Gill_Vwnd(IX,IY) = AllFrc_SrfFlux.AllFrc_Gill_Vwnd(IX,IY) + &
          ( EarthBeta * GrdInfo.VReY(IY) ) * (-0.25d0) * &
          ( SrfFlux.Gill_Uwnd(IX+0,IY-1) + SrfFlux.Gill_Uwnd(IX+0,IY+0) + &
            SrfFlux.Gill_Uwnd(IX+1,IY-1) + SrfFlux.Gill_Uwnd(IX+1,IY+0) );
          
          !!! Vwnd Damping
          AllFrc_SrfFlux.AllFrc_Gill_Vwnd(IX,IY) = AllFrc_SrfFlux.AllFrc_Gill_Vwnd(IX,IY) + &
          (-1.0d0) * (Gill_UVwnd_Damp) * SrfFlux.Gill_Vwnd(IX,IY)
    
         END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IX,IY)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO IX = 1,NX
       DO IY = 1,NY
       
           !!! C^2 * (Ux)
           AllFrc_SrfFlux.AllFrc_Gill_Pres(IX+1,IY+1) = &
           (-1.0d0/DX) * (AtmEigC_Pow2) * ( SrfFlux.Gill_Uwnd(IX+1,IY) - SrfFlux.Gill_Uwnd(IX,IY) );
           
           !!! C^2 *(Vy)
           AllFrc_SrfFlux.AllFrc_Gill_Pres(IX+1,IY+1) = AllFrc_SrfFlux.AllFrc_Gill_Pres(IX+1,IY+1) + &
           (-1.0d0/DY) * (AtmEigC_Pow2) * (SrfFlux.Gill_Vwnd(IX,IY+1) - SrfFlux.Gill_Vwnd(IX,IY) );
           
           !!! Pres Damping
           AllFrc_SrfFlux.AllFrc_Gill_Pres(IX+1,IY+1) = AllFrc_SrfFlux.AllFrc_Gill_Pres(IX+1,IY+1) + &
           (-1.0d0) * (Gill_Pres_Damp) * SrfFlux.Gill_Pres(IX+1,IY+1)
           
           !!!! SSTA Heating
           AllFrc_SrfFlux.AllFrc_Gill_Pres(IX+1,IY+1) = AllFrc_SrfFlux.AllFrc_Gill_Pres(IX+1,IY+1) + &
           (-1.0d0) * (Gill_AtmBLCoeff) * (OcnStat.SST(IX+1,IY+1)) * &
                       exp( (Clim_OcnStat.SSTBar(IX+1,IY+1,CMonth)-30.0d0-273.15d0)/16.7d0 );
       
       END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

END SUBROUTINE Get_AllFrc_SurfaceFlux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Add_AllFrc_SurfaceFlux &
           (GrdInfo,Clim_OcnStat, TimTick, Old_SrfFlux, New_SrfFlux, AllFrc_SrfFlux)

    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(AllFrc_SurfaceFlux) :: AllFrc_SrfFlux
    TYPE(SurfaceFlux) :: Old_SrfFlux, New_SrfFlux
    
    integer*8 :: NX,NY,NZ,NMode
    integer*8 :: IX,IY,IZ,IMode
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    NZ = GrdInfo.NZ
    NMode = GrdInfo.NMode
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IX,IY)
    DO IX = 1,NX+1
       DO IY = 1,NY
       
       New_SrfFlux.Gill_Uwnd(IX,IY) = Old_SrfFlux.Gill_Uwnd(IX,IY) + &
       (1.0d0/dble(AtmBLSplit)) * (GrdInfo.ReDT) * AllFrc_SrfFlux.AllFrc_Gill_Uwnd(IX,IY);
       
       END DO
    END DO
    !$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IX,IY)
    DO IX = 1,NX
       DO IY = 1,NY+1
       
       New_SrfFlux.Gill_Vwnd(IX,IY) = &
       Old_SrfFlux.Gill_Vwnd(IX,IY) + &
       (1.0d0/dble(AtmBLSplit)) * (GrdInfo.ReDT) * AllFrc_SrfFlux.AllFrc_Gill_Vwnd(IX,IY)
       
       END DO
    END DO
    !$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IX,IY)
    DO IX = 1,NX+2
       DO IY = 1,NY+2
       
       New_SrfFlux.Gill_Pres(IX,IY) = &
       Old_SrfFlux.Gill_Pres(IX,IY) + &
       (1.0d0/dble(AtmBLSplit)) * (GrdInfo.ReDT) * AllFrc_SrfFlux.AllFrc_Gill_Pres(IX,IY)
       
       END DO
    END DO
    !$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! boundary condition
    New_SrfFlux.Gill_Uwnd(1,1:NY) = 0.0d0;
    New_SrfFlux.Gill_Uwnd(NX+1,1:NY) = 0.0d0;
    
    New_SrfFlux.Gill_Vwnd(1:NX,1) = 0.0d0;
    New_SrfFlux.Gill_Vwnd(1:NX,NY+1) = 0.0d0;
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE Add_AllFrc_SurfaceFlux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OneStep_Forward_Gill &
           (GrdInfo,Clim_OcnStat, TimTick, SrfFlux, Tmp_SrfFlux, AllFrc_SrfFlux,OcnStat)

    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(AllFrc_SurfaceFlux) :: AllFrc_SrfFlux
    TYPE(SurfaceFlux) :: SrfFlux, Tmp_SrfFlux
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!! Matsuno Scheme 
    !!!! F(Un)
    call Get_AllFrc_SurfaceFlux &
         (GrdInfo,Clim_OcnStat, TimTick, SrfFlux, AllFrc_SrfFlux, OcnStat)
    
    !!! U*  = Un + DT * F(Un)
    call Add_AllFrc_SurfaceFlux &
        (GrdInfo,Clim_OcnStat, TimTick, SrfFlux, Tmp_SrfFlux, AllFrc_SrfFlux)
    
    !!!! F(U*) 
    call Get_AllFrc_SurfaceFlux &
        (GrdInfo,Clim_OcnStat, TimTick, Tmp_SrfFlux, AllFrc_SrfFlux, OcnStat)
    
    !!! Un+1 = Un + DT * F(U*)
    call Add_AllFrc_SurfaceFlux &
        (GrdInfo,Clim_OcnStat, TimTick, SrfFlux, SrfFlux, AllFrc_SrfFlux)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!! Euler Scheme
    !!!!! F(Un)
    !call Get_AllFrc_SurfaceFlux &
    !     (GrdInfo,Clim_OcnStat, TimTick, SrfFlux, AllFrc_SrfFlux, OcnStat)
    !
    !!!! Un+1 = Un + DT * F(U*)
    !call Add_AllFrc_SurfaceFlux &
    !    (GrdInfo,Clim_OcnStat, TimTick, SrfFlux, SrfFlux, AllFrc_SrfFlux)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE OneStep_Forward_Gill
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OneStep_Forward_Gill_Lindzen &
           (GrdInfo,Clim_OcnStat, TimTick, SrfFlux, Tmp_SrfFlux, AllFrc_SrfFlux,OcnStat)

    !!! Linearized Lindzen-Nigam Response 

    !!! (Ra)*U - (Beta*Y*V) = (dR)*(Tx) = Fx
    !!! (Ra)*V + (Beta*Y*U) = (dR)*(Ty) = Fy

    !!! U  = (Ra*Fx + Beta*Y*Fy)/(Ra^2 + (Beta*Y)^2);
    !!! V  = (Ra*Fy - Beta*Y*Fx)/(Ra^2 + (Beta*Y)^2);
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(AllFrc_SurfaceFlux) :: AllFrc_SrfFlux
    TYPE(SurfaceFlux) :: SrfFlux, Tmp_SrfFlux
    
    real*8,allocatable :: pSSTA_pX(:,:),pSSTA_pY(:,:)

    integer*8 :: NX,NY
    integer*8 :: IX,IY
    real*8 :: DX,DY
    
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    
    DX = GrdInfo.ReDX
    DY = GrdInfo.ReDY
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! ArakaWa C Grid, Equi-Distant Nodes 
    !!!
    !!!                                        U(IX-1,IY+0)
    !!!              P(IX+0,IY+0) V(IX-1,IY+0) P(IX+0,IY+1) V(IX-1,IY+1) P(IX+0,IY+2)
    !!!              U(IX+0,IY-1)              U(IX+0,IY+0)              U(IX+0,IY+1)
    !!! V(IX+0,IY-1) P(IX+1,IY+0) V(IX+0,IY+0) P(IX+1,IY+1) V(IX+0,IY+1) P(IX+1,IY+2) V(IX+0,IY+2)
    !!!              U(IX+1,IY-1)              U(IX+1,IY+0)              U(IX+1,IY+1)
    !!!              P(IX+2,IY+0) V(IX+1,IY+0) P(IX+2,IY+1) V(IX+1,IY+1) P(IX+2,IY+2)
    !!!                                        U(IX+2,IY+0)                      
    !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IX,IY)
    DO IX = 1,NX+1
       DO IY = 1,NY
       
        SrfFlux.Gill_Uwnd(IX,IY) = & 
        Gill_UVwnd_Damp * (1.0d0/DX) * (OcnStat.SST(IX+1,IY+1) - OcnStat.SST(IX+0,IY+1) );
       
        
        SrfFlux.Gill_Uwnd(IX,IY) = SrfFlux.Gill_Uwnd(IX,IY) + &
        (+1.0d0/(2.0d0*DY) ) *( (0.5d0)*(OcnStat.SST(IX+0,IY+2) + OcnStat.SST(IX+1,IY+2)) - &
                                (0.5d0)*(OcnStat.SST(IX+0,IY+0) + OcnStat.SST(IX+1,IY+0)) );
        
        SrfFlux.Gill_Uwnd(IX,IY) = SrfFlux.Gill_Uwnd(IX,IY)/ &
        ( Gill_UVwnd_Damp**(2.0d0) + (EarthBeta * GrdInfo.UReY(IY))**(2.0d0) );
       
        
       END DO
    END DO
    !$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IX,IY)
    DO IX = 1,NX
       DO IY = 1,NY+1
       
        SrfFlux.Gill_Vwnd(IX,IY) = &
        Gill_UVwnd_Damp * (1.0d0/DY) * (OcnStat.SST(IX+1,IY+1) - OcnStat.SST(IX+1,IY+0) );
        
        SrfFlux.Gill_Vwnd(IX,IY) = SrfFlux.Gill_Vwnd(IX,IY) + &
        (-1.0d0/(2.0d0*DX) ) * ( (0.5d0)*(OcnStat.SST(IX+2,IY+1) + OcnStat.SST(IX+2,IY+0)) - &
                                 (0.5d0)*(OcnStat.SST(IX+0,IY+1) + OcnStat.SST(IX+0,IY+0)) );
          
        SrfFlux.Gill_Vwnd(IX,IY) = SrfFlux.Gill_Vwnd(IX,IY)/ &
        ( Gill_UVwnd_Damp**(2.0d0) + (EarthBeta * GrdInfo.VReY(IY))**(2.0d0) );
       
       
       END DO
    END DO
    !$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   

END SUBROUTINE OneStep_Forward_Gill_Lindzen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE MO_ATMBLSOLVER
