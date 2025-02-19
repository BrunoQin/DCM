MODULE MO_OCEANSOLVER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Numerical Solver for Ocean 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE MO_NCTOOLS
USE MO_MATHTOOLS
USE MO_TYPEDEF
USE MO_TAUTSUBMODEL

IMPLICIT NONE

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Get_AllFrc_OceanStat_OneMode &
           (GrdInfo,Clim_OcnStat, TimTick, OcnStat, AllFrc_OcnStat, SrfFlux,IM)

!!!! get one-baroclinic mode force, Equivalent Barotropic Solver
!!!! Shallow Water Solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat, Tmp_OcnStat 
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(AllFrc_OceanStat) :: AllFrc_OcnStat
    TYPE(SurfaceFlux) :: SrfFlux
    
    integer*8 :: NMode
    integer*8 :: IM !!! indicate the index of baroclinic mode
    
    integer*8 :: NX,NY,NZ
    integer*8 :: IX,IY,IZ
    
    integer*8 :: CMonth !!! the calendar month
    
    real*8 :: DX,DY
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    NZ = GrdInfo.NZ
    NMode = GrdInfo.NMode
    
    DX = GrdInfo.ReDX
    DY = GrdInfo.ReDY
    
    CMonth = TimTick.CMonth
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IX,IY)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! F(Ucur)
    DO IX = 2,NX
       DO IY = 1,NY
       
          !!! pressure gradient force
          AllFrc_OcnStat.AllFrc_Co_Ucur(IX,IY,IM) = &
          (-1.d0/DX) * ( OcnStat.Co_Pres(IX+1,IY+1,IM) - OcnStat.Co_Pres(IX,IY+1,IM) );
          
          !!!! Coriolis force
          AllFrc_OcnStat.AllFrc_Co_Ucur(IX,IY,IM) = AllFrc_OcnStat.AllFrc_Co_Ucur(IX,IY,IM) + &
          ( EarthBeta * GrdInfo.UReY(IY) ) * (0.25d0) * &
          ( OcnStat.Co_Vcur(IX-1,IY+0,IM) + OcnStat.Co_Vcur(IX-1,IY+1,IM) + &
            OcnStat.Co_Vcur(IX+0,IY+0,IM) + OcnStat.Co_Vcur(IX+0,IY+1,IM) );
          
          !!! horizontal eddy dissipation for Ucur
           AllFrc_OcnStat.AllFrc_Co_Ucur(IX,IY,IM) = AllFrc_OcnStat.AllFrc_Co_Ucur(IX,IY,IM) + &
          ( XY_Diff_Co_UVcur ) * (1.0d0/(DX*DX)) * ( OcnStat.Co_Ucur(IX-1,IY,IM) + &
                                                     OcnStat.Co_Ucur(IX+1,IY,IM) + &
                                                     (-2.0d0) * OcnStat.Co_Ucur(IX,IY,IM) )
          
          AllFrc_OcnStat.AllFrc_Co_Ucur(IX,IY,IM) = AllFrc_OcnStat.AllFrc_Co_Ucur(IX,IY,IM) + &
          ( XY_Diff_Co_UVcur ) * (1.0d0/(DY*DY)) * ( OcnStat.Co_Ucur(IX,safe(IY-1,NY),IM) + &
                                                     OcnStat.Co_Ucur(IX,safe(IY+1,NY),IM) + &
                                                     (-2.0d0) * OcnStat.Co_Ucur(IX,IY,IM) )
          
          !!! vertical dissipation (equivalent damping) for Ucur
          AllFrc_OcnStat.AllFrc_Co_Ucur(IX,IY,IM) = AllFrc_OcnStat.AllFrc_Co_Ucur(IX,IY,IM) + &
          ( Z_Diff_Co_UVcur/Clim_OcnStat.EigC_Pow2(IX+1,IM)) *(-1.0d0) * OcnStat.Co_Ucur(IX,IY,IM)
          
          
          !!!! X-direction wind stress
          AllFrc_OcnStat.AllFrc_Co_Ucur(IX,IY,IM) = AllFrc_OcnStat.AllFrc_Co_Ucur(IX,IY,IM) + &
          ( Clim_OcnStat.Co_Wstr(IX+1,IM)/RhoWat ) * (Wstr_Proj_Amp) * SrfFlux.TAUX(IX,IY);          
          
       
       END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IX,IY)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! F(Vcur)
    DO IX = 1,NX
       DO IY = 2,NY
       
          !!! pressure gradient force
          AllFrc_OcnStat.AllFrc_Co_Vcur(IX,IY,IM) = &
          (-1.d0/DY) * ( OcnStat.Co_Pres(IX+1,IY+1,IM) - OcnStat.Co_Pres(IX+1,IY,IM) );
          
          !!!!! Coriolis force
          AllFrc_OcnStat.AllFrc_Co_Vcur(IX,IY,IM) = AllFrc_OcnStat.AllFrc_Co_Vcur(IX,IY,IM) + &
          ( EarthBeta * GrdInfo.VReY(IY) ) * (-0.25d0) * &
          ( OcnStat.Co_Ucur(IX+0,IY-1,IM) + OcnStat.Co_Ucur(IX+0,IY+0,IM) + &
            OcnStat.Co_Ucur(IX+1,IY-1,IM) + OcnStat.Co_Ucur(IX+1,IY+0,IM) );
          
          !!! horizontal eddy dissipation for Vcur
          AllFrc_OcnStat.AllFrc_Co_Vcur(IX,IY,IM) = AllFrc_OcnStat.AllFrc_Co_Vcur(IX,IY,IM) + &
          ( XY_Diff_Co_UVcur ) * (1.0d0/(DX*DX)) * ( OcnStat.Co_Vcur(safe(IX-1,NX),IY,IM) + &
                                                     OcnStat.Co_Vcur(safe(IX+1,NX),IY,IM) + &
                                                     (-2.0d0) * OcnStat.Co_Vcur(IX,IY,IM) )
          
          AllFrc_OcnStat.AllFrc_Co_Vcur(IX,IY,IM) = AllFrc_OcnStat.AllFrc_Co_Vcur(IX,IY,IM) + &
          ( XY_Diff_Co_UVcur ) * (1.0d0/(DY*DY)) * ( OcnStat.Co_Vcur(IX,IY-1,IM) + &
                                                     OcnStat.Co_Vcur(IX,IY+1,IM) + &
                                                     (-2.0d0) * OcnStat.Co_Vcur(IX,IY,IM) )
          
          !!! vertical dissipation (equivalent damping) for Vcur
          AllFrc_OcnStat.AllFrc_Co_Vcur(IX,IY,IM) = AllFrc_OcnStat.AllFrc_Co_Vcur(IX,IY,IM) + &
          ( Z_Diff_Co_UVcur/Clim_OcnStat.EigC_Pow2(IX+1,IM)) *(-1.0d0) * OcnStat.Co_Vcur(IX,IY,IM)
          
          !!!! Y-direction wind stress
          AllFrc_OcnStat.AllFrc_Co_Vcur(IX,IY,IM) = AllFrc_OcnStat.AllFrc_Co_Vcur(IX,IY,IM) + &
          ( Clim_OcnStat.Co_Wstr(IX+1,IM) / RhoWat ) * (Wstr_Proj_Amp) * SrfFlux.TAUY(IX,IY);
       
       END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IX,IY)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! F(Pres)
    DO IX = 1,NX
        DO IY = 1,NY
        
           !!! C^2 * (Ux)
           AllFrc_OcnStat.AllFrc_Co_Pres(IX+1,IY+1,IM) = &
           (-1.0d0/DX) * (Clim_OcnStat.EigC_Pow2(IX+1,IM))*  &
                         ( OcnStat.Co_Ucur(IX+1,IY,IM) - OcnStat.Co_Ucur(IX,IY,IM) );
           
           !!! C^2 *(Vy)
           AllFrc_OcnStat.AllFrc_Co_Pres(IX+1,IY+1,IM) = AllFrc_OcnStat.AllFrc_Co_Pres(IX+1,IY+1,IM) + &
           (-1.0d0/DY) * (Clim_OcnStat.EigC_Pow2(IX+1,IM))*  &
                         ( OcnStat.Co_Vcur(IX,IY+1,IM) - OcnStat.Co_Vcur(IX,IY,IM) );
           
           
           !!!! horizontal dissipation for Pres
           AllFrc_OcnStat.AllFrc_Co_Pres(IX+1,IY+1,IM) = AllFrc_OcnStat.AllFrc_Co_Pres(IX+1,IY+1,IM) + &
           ( XY_Diff_Co_Pres ) * (1.0d0/(DX*DX)) * ( OcnStat.Co_Pres(IX,IY+1,IM) + &
                                                     (-2.0d0) * OcnStat.Co_Pres(IX+1,IY+1,IM) + &
                                                     OcnStat.Co_Pres(IX+2,IY+1,IM) )
           
           AllFrc_OcnStat.AllFrc_Co_Pres(IX+1,IY+1,IM) = AllFrc_OcnStat.AllFrc_Co_Pres(IX+1,IY+1,IM) + &
           ( XY_Diff_Co_Pres ) * (1.0d0/(DY*DY)) * ( OcnStat.Co_Pres(IX+1,IY,IM) + &
                                                     (-2.0d0) * OcnStat.Co_Pres(IX+1,IY+1,IM) + &
                                                     OcnStat.Co_Pres(IX+1,IY+2,IM) )
           
           !!! vertical dissipation (equivalent damping) for pres
          AllFrc_OcnStat.AllFrc_Co_Pres(IX+1,IY+1,IM) = AllFrc_OcnStat.AllFrc_Co_Pres(IX+1,IY+1,IM) + &
          ( Z_Diff_Co_Pres/Clim_OcnStat.EigC_Pow2(IX+1,IM)) *(-1.0d0) * OcnStat.Co_Pres(IX+1,IY+1,IM)
              
        END DO
     END DO
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !$omp end parallel do
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

END SUBROUTINE Get_AllFrc_OceanStat_OneMode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Get_AllFrc_OceanStat_MulMode &
           (GrdInfo,Clim_OcnStat, TimTick, OcnStat, AllFrc_OcnStat,SrfFlux)

!!!! get multi-baroclinic mode force

    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat, Tmp_OcnStat 
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(AllFrc_OceanStat) :: AllFrc_OcnStat
    TYPE(SurfaceFlux) :: SrfFlux
    
    integer*8 :: IM
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO IM = 1, GrdInfo.NMode
    
        call Get_AllFrc_OceanStat_OneMode &
             (GrdInfo,Clim_OcnStat, TimTick, OcnStat, AllFrc_OcnStat, SrfFlux, IM)
    
    END DO !!!! IMode
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE Get_AllFrc_OceanStat_MulMode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Add_AllFrc_OceanStat_MulMode &
           (GrdInfo,Clim_OcnStat, TimTick, Old_OcnStat, New_OcnStat, AllFrc_OcnStat)

     !!! New_OcnStat = Old_OcnStat + (GrdInfo.ReDT) * AllFrc_OcnStat (Euler scheme)
     !!! Mask

     TYPE(GridInfo) :: GrdInfo
     TYPE(Clim_OceanStat) :: Clim_OcnStat
     TYPE(TimeTick) :: TimTick
     
     TYPE(OceanStat) :: Old_OcnStat, New_OcnStat 
    
     TYPE(AllFrc_OceanStat) :: AllFrc_OcnStat
     
     integer*8 :: NX,NY,NZ,NMode
     integer*8 :: IX,IY,IZ,IMode
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     NX = GrdInfo.NX
     NY = GrdInfo.NY
     NZ = GrdInfo.NZ
     NMode = GrdInfo.NMode
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !$omp parallel do schedule (static) &
     !$omp default(shared) &
     !$omp private(IMODE,IX,IY)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!! Add AllFrc_Co_Pres & Mask
     DO IMode = 1,NMode
        DO IX = 1,NX+2
           DO IY = 1,NY+2
           
            New_OcnStat.Co_Pres(IX,IY,IMode) = &
            Old_OcnStat.Co_Pres(IX,IY,IMode) + (GrdInfo.ReDT) * AllFrc_OcnStat.AllFrc_Co_Pres(IX,IY,IMode)
            
            New_OcnStat.Co_Pres(IX,IY,IMode) = New_OcnStat.Co_Pres(IX,IY,IMode) * GrdInfo.PMask(IX,IY)

           END DO
        END DO
     END DO
     !$omp end parallel do
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !$omp parallel do schedule (static) &
     !$omp default(shared) &
     !$omp private(IMODE,IX,IY)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!! Add AllFrc_Co_Ucur & Mask
     DO IMode = 1,NMode
        DO IX = 1,NX+1
           DO IY = 1,NY
           
            New_OcnStat.Co_Ucur(IX,IY,IMode) = &
            Old_OcnStat.Co_Ucur(IX,IY,IMode) + (GrdInfo.ReDT) * AllFrc_OcnStat.AllFrc_Co_Ucur(IX,IY,IMode)
            
            New_OcnStat.Co_Ucur(IX,IY,IMode) = New_OcnStat.Co_Ucur(IX,IY,IMode) * GrdInfo.UMask(IX,IY)

           END DO
        END DO
     END DO
     !$omp end parallel do
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !$omp parallel do schedule (static) &
     !$omp default(shared) &
     !$omp private(IMODE,IX,IY)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!! Add AllFrc_Co_Vcur & Mask
     DO IMode = 1,NMode
        DO IX = 1,NX
           DO IY = 1,NY+1
           
            New_OcnStat.Co_Vcur(IX,IY,IMode) = &
            Old_OcnStat.Co_Vcur(IX,IY,IMode) + (GrdInfo.ReDT) * AllFrc_OcnStat.AllFrc_Co_Vcur(IX,IY,IMode)
            
            New_OcnStat.Co_Vcur(IX,IY,IMode) = New_OcnStat.Co_Vcur(IX,IY,IMode) * GrdInfo.VMask(IX,IY)

           END DO
        END DO
     END DO
     !$omp end parallel do
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     
     !!! Sponge Boundary Condition at North/South boundary
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !$omp parallel do schedule (static) &
     !$omp default(shared) &
     !$omp private(IMODE,IY)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     DO IMode = 1,NMode
        DO IY = 1, GrdInfo.Sponge_Num
        
        !!! Ucur
        New_OcnStat.Co_Ucur(1:NX+1,IY,IMode) = &
        (1.0d0 - GrdInfo.Sponge_Alp(IY) ) * New_OcnStat.Co_Ucur(1:NX+1,IY,IMode)
        
        New_OcnStat.Co_Ucur(1:NX+1,NY-(IY-1),IMode) = &
        (1.0d0 - GrdInfo.Sponge_Alp(IY) ) * New_OcnStat.Co_Ucur(1:NX+1,NY-(IY-1),IMode)
     
        !!! Vcur
        New_OcnStat.Co_Vcur(1:NX,IY,IMode) = &
        (1.0d0 - GrdInfo.Sponge_Alp(IY) ) * New_OcnStat.Co_Vcur(1:NX,IY,IMode)
        
        New_OcnStat.Co_Vcur(1:NX,NY+1-(IY-1),IMode) = &
        (1.0d0 - GrdInfo.Sponge_Alp(IY) ) * New_OcnStat.Co_Vcur(1:NX,NY+1-(IY-1),IMode)
        
        !!! Pres
        New_OcnStat.Co_Pres(1:NX+2,IY,IMode) = &
        (1.0d0 - GrdInfo.Sponge_Alp(IY) ) * New_OcnStat.Co_Pres(1:NX+2,IY,IMode)
        
        New_OcnStat.Co_Pres(1:NX+2,NY+2-(IY-1),IMode) = &
        (1.0d0 - GrdInfo.Sponge_Alp(IY) ) * New_OcnStat.Co_Pres(1:NX+2,NY+2-(IY-1),IMode)
     
        END DO
     END DO
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !$omp end parallel do
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE Add_AllFrc_OceanStat_MulMode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Get_AllFrc_OceanStat_SSTA &
           (GrdInfo,Clim_OcnStat, TimTick, OcnStat, AllFrc_OcnStat,SrfFlux)

    !!!! get multi-baroclinic mode force

    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat, Tmp_OcnStat 
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(AllFrc_OceanStat) :: AllFrc_OcnStat
    TYPE(SurfaceFlux) :: SrfFlux
    
    integer*8 :: NX,NY,NZ,NMode
    integer*8 :: IX,IY,IZ,IMode
    integer*8 :: CMonth !!! the calendar month
    integer*8 :: YrDaySum 
    
    real*8 :: DX,DY,DT
    
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    NZ = GrdInfo.NZ
    NMode = GrdInfo.NMode
    
    DX = GrdInfo.ReDX
    DY = GrdInfo.ReDY
    DT = GrdInfo.ReDT
    
    CMonth = TimTick.CMonth
    YrDaySum = TimTick.YrDaySum
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IX,IY)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO IX = 1,NX
       DO IY = 1,NY
       
          
          !!! UATBX (Zonal Advective FeedBack) 
          AllFrc_OcnStat.UATBX(IX+1,IY+1) = (UATBX_Amp) * &
          UpWindXYAdv(UL = OcnStat.Usrf(IX,IY), &
                      UR = OcnStat.Usrf(IX+1,IY), &
                      TL = Clim_OcnStat.SSTBar(IX,IY+1,CMonth), &
                      TM = Clim_OcnStat.SSTBar(IX+1,IY+1,CMonth), &
                      TR = Clim_OcnStat.SSTBar(IX+2,IY+1,CMonth),DXY = DX)
          
          !!!! UBTAX
          AllFrc_OcnStat.UBTAX(IX+1,IY+1) = (UBTAX_Amp) * &
          UpWindXYAdv(UL = Clim_OcnStat.UsrfBar(IX,IY,CMonth), &
                      UR = Clim_OcnStat.UsrfBar(IX+1,IY,CMonth), &
                      TL = OcnStat.SST(IX,IY+1), &
                      TM = OcnStat.SST(IX+1,IY+1), &
                      TR = OcnStat.SST(IX+2,IY+1),DXY = DX)
          
          
          !!! UATAX
          AllFrc_OcnStat.UATAX(IX+1,IY+1) = (UATAX_Amp) * &
          UpWindXYAdv(UL = OcnStat.Usrf(IX,IY), &
                      UR = OcnStat.Usrf(IX+1,IY), &
                      TL = OcnStat.SST(IX,IY+1), &
                      TM = OcnStat.SST(IX+1,IY+1), &
                      TR = OcnStat.SST(IX+2,IY+1),DXY = DX)
          
          
          !!! VATBY
          AllFrc_OcnStat.VATBY(IX+1,IY+1) = (VATBY_Amp) * &
          UpWindXYAdv(UL = OcnStat.Vsrf(IX,IY), &
                      UR = OcnStat.Vsrf(IX,IY+1), &
                      TL = Clim_OcnStat.SSTBar(IX+1,IY,CMonth), &
                      TM = Clim_OcnStat.SSTBar(IX+1,IY+1,CMonth), &
                      TR = Clim_OcnStat.SSTBar(IX+1,IY+2,CMonth),DXY = DY)
          
          !!! VBTAY
          AllFrc_OcnStat.VBTAY(IX+1,IY+1) = (VBTAY_Amp) * &
          UpWindXYAdv(UL = Clim_OcnStat.VsrfBar(IX,IY,CMonth), &
                      UR = Clim_OcnStat.VsrfBar(IX,IY+1,CMonth), &
                      TL = OcnStat.SST(IX+1,IY), &
                      TM = OcnStat.SST(IX+1,IY+1), &
                      TR = OcnStat.SST(IX+1,IY+2),DXY = DY)
          
          !!! VATAY
          AllFrc_OcnStat.VATAY(IX+1,IY+1) =  (VATAY_Amp) * &
          UpWindXYAdv(UL = OcnStat.Vsrf(IX,IY), &
                      UR = OcnStat.Vsrf(IX,IY+1), &
                      TL = OcnStat.SST(IX+1,IY), &
                      TM = OcnStat.SST(IX+1,IY+1), &
                      TR = OcnStat.SST(IX+1,IY+2),DXY = DY)
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!! - [M(WA+WB)(TA+TB)z - M(WB)(TB)z] = &
          !!! - (M(WA+WB)-M(WB))*(TB)z - (M(WA+WB)-M(WB))*(TA)z - M(WB)*(TA)z
          
          !!! Ekman FeedBack       = WATBZ = - (M(WA+WB)-M(WB))*(TB)z
          !!! Thermocline FeedBack = WBTAZ = - M(WB)*(TA)z
          !!! Nonlinear Term       = WATAZ = - (M(WA+WB)-M(WB))*(TA)z
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          !!!! WATBZ (Ekman FeedBack)
          !!! let mixed layer depth ~= 50.0m or 75.0m 
          AllFrc_OcnStat.WATBZ(IX+1,IY+1) = (WATBZ_Amp) * (-1.0d0) * &
          ( psgn( Clim_OcnStat.WsrfBar(IX+1,IY+1,CMonth) + OcnStat.Wsrf(IX+1,IY+1) )  - &
            psgn( Clim_OcnStat.WsrfBar(IX+1,IY+1,CMonth) ) ) * &
          ( (Clim_OcnStat.SSTBar(IX+1,IY+1,CMonth)-Clim_OcnStat.TSUBBar(IX+1,IY+1,CMonth))/ 50.0d0 )
          
          !!! WBTAZ (Thermocline FeedBack)
          AllFrc_OcnStat.WBTAZ(IX+1,IY+1) = (WBTAZ_Amp) * (-1.0d0) * (Clim_OcnStat.THFK_Mask(IX+1,IY+1)) * &
          psgn( Clim_OcnStat.WsrfBar(IX+1,IY+1,CMonth) ) *  &
          ( (OcnStat.SST(IX+1,IY+1) - OcnStat.TSUB(IX+1,IY+1))/(50.0d0) )
          
          !!! WATAZ (nonlinear vertical heating)
           AllFrc_OcnStat.WATAZ(IX+1,IY+1) = (WATAZ_Amp) * (-1.0d0) * (Clim_OcnStat.THFK_Mask(IX+1,IY+1)) * &
          ( psgn( Clim_OcnStat.WsrfBar(IX+1,IY+1,CMonth) + OcnStat.Wsrf(IX+1,IY+1) ) - &
            psgn( Clim_OcnStat.WsrfBar(IX+1,IY+1,CMonth) ) ) * &
          ( (OcnStat.SST(IX+1,IY+1) - OcnStat.TSUB(IX+1,IY+1))/(50.0d0) )

          !!! SSTA_RESQ (SSTA Diffusion)
          AllFrc_OcnStat.SSTA_RESQ(IX+1,IY+1) =  &
          (DiffX_SSTA) * &
          ( OcnStat.SST(IX,IY+1) - 2.0d0*OcnStat.SST(IX+1,IY+1) + OcnStat.SST(IX+2,IY+1) )/DX**(2.0d0) + &
          (DiffY_SSTA) * &
          ( OcnStat.SST(IX+1,IY) - 2.0d0*OcnStat.SST(IX+1,IY+1) + OcnStat.SST(IX+1,IY+2) )/DY**(2.0d0)

          
          !!! SSTA_RESQ (seasonal SSTA Damping)
          AllFrc_OcnStat.SSTA_RESQ(IX+1,IY+1) =  AllFrc_OcnStat.SSTA_RESQ(IX+1,IY+1) +  &
          (-Damp_SSTA) * (1.0d0 + (PhaseLocking_OnOff) * sin( 2.0d0 * (pi/365.0d0) * ( dble(YrDaySum) + Phase_Difference) ) ) * &
          (OcnStat.SST(IX+1,IY+1))**(SSTA_DampPower) 
          
          
          !!!! SSTA_RESQ (SSTA Vertical Diffusion)
          AllFrc_OcnStat.SSTA_RESQ(IX+1,IY+1) = AllFrc_OcnStat.SSTA_RESQ(IX+1,IY+1) +  &
          ( DiffZ_SSTA ) * ( OcnStat.TSUB (IX+1,IY+1) - OcnStat.SST(IX+1,IY+1) )/( (50.0d0) * (100.0d0) );
          
          
          !!! Compute ZAFK, THFK, RESFK by the way
          AllFrc_OcnStat.ZAFK(IX+1,IY+1) = AllFrc_OcnStat.UATBX(IX+1,IY+1);
          
          AllFrc_OcnStat.THFK(IX+1,IY+1) = AllFrc_OcnStat.WBTAZ(IX+1,IY+1);
          
          AllFrc_OcnStat.RESFK(IX+1,IY+1) = AllFrc_OcnStat.UBTAX(IX+1,IY+1) + &
                                            AllFrc_OcnStat.UATAX(IX+1,IY+1) + &
                                            AllFrc_OcnStat.VATBY(IX+1,IY+1) + &
                                            AllFrc_OcnStat.VBTAY(IX+1,IY+1) + &
                                            AllFrc_OcnStat.VATAY(IX+1,IY+1) + &
                                            AllFrc_OcnStat.WATBZ(IX+1,IY+1) + &
                                            AllFrc_OcnStat.WATAZ(IX+1,IY+1) + &
                                            AllFrc_OcnStat.SSTA_RESQ(IX+1,IY+1)
       
       END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!! Nudging of SSTA 
    IF ( ( SSTA_Nudging_OnOff == 1 ).AND.(Eq_Direct_Assign == 0) ) THEN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IX,IY)
    DO IX = 1,NX+2
       DO IY = 1,NY+2
    
       IF( (GrdInfo.PLon(IX) >= Nudging_WLon ).AND. ( GrdInfo.PLon(IX) <= Nudging_ELon ) )THEN
          IF( (GrdInfo.PLat(IY) >= Nudging_SLat ).AND. ( GrdInfo.PLat(IY) <= Nudging_NLat ) )THEN
           
           !!! NUDGE TERM
           AllFrc_OcnStat.NUDGE(IX,IY) = ( (1.0d0)/(SSTA_Nudging_Scale * Day_DT) ) * &
                                         ( OcnStat.Nudging_SST(IX,IY) - OcnStat.SST(IX,IY) )

           !!! Change Merdional Nudging Amplitude
           !!! AllFrc_OcnStat.NUDGE(IX,IY) = &
           !!! ( 1.0d0 + min( abs(GrdInfo.PLat(IY)),15.0d0)/2.5d0 ) * AllFrc_OcnStat.NUDGE(IX,IY) 

           !!! Influence Radius
           !!! AllFrc_OcnStat.NUDGE(IX,IY) = &
           !!! exp( -(GrdInfo.PLat(IY)/SSTA_Nudging_Radius)**(2.0d0) ) * AllFrc_OcnStat.NUDGE(IX,IY);


           !!! ANOTHER NUDGING 
           !AllFrc_OcnStat.NUDGE(IX,IY) = (WATAZ_Amp) * (-1.0d0) * (Clim_OcnStat.THFK_Mask(IX,IY)) * &
           !( psgn( Clim_OcnStat.WsrfBar(IX,IY,CMonth) + OcnStat.Wsrf(IX,IY) ) )* &
           !( (OcnStat.SST(IX,IY) - OcnStat.Nudging_SST(IX,IY))/(50.0d0) ) 
          
           
           !!! RESQ TERM = RESQ TERM + NUDGE TERM
           AllFrc_OcnStat.SSTA_RESQ(IX,IY) = AllFrc_OcnStat.SSTA_RESQ(IX,IY) + &
                                             AllFrc_OcnStat.NUDGE(IX,IY)
           

          END IF
       END IF 
    
       END DO
    END DO
    !$omp end parallel do
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
END SUBROUTINE Get_AllFrc_OceanStat_SSTA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Add_AllFrc_OceanStat_SSTA &
           (GrdInfo,Clim_OcnStat, TimTick, Old_OcnStat, New_OcnStat, AllFrc_OcnStat)

     !!! New_OcnStat = Old_OcnStat + (GrdInfo.ReDT) * AllFrc_OcnStat (Euler scheme)
     !!! Mask

     TYPE(GridInfo) :: GrdInfo
     TYPE(Clim_OceanStat) :: Clim_OcnStat
     TYPE(TimeTick) :: TimTick
     
     TYPE(OceanStat) :: Old_OcnStat, New_OcnStat 
    
     TYPE(AllFrc_OceanStat) :: AllFrc_OcnStat
     
     integer*8 :: NX,NY,NZ,NMode
     integer*8 :: IX,IY,IZ,IMode
     
     integer*8 :: CMonth !!! the calendar month
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     NX = GrdInfo.NX
     NY = GrdInfo.NY
     NZ = GrdInfo.NZ
     NMode = GrdInfo.NMode
     
     CMonth = TimTick.CMonth
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !$omp parallel do schedule (static) &
     !$omp default(shared) &
     !$omp private(IX,IY)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     DO IX = 1,NX+2
           DO IY = 1,NY+2
           
            !!! UATBX
            New_OcnStat.SST(IX,IY) = &
            Old_OcnStat.SST(IX,IY) + (GrdInfo.ReDT) * AllFrc_OcnStat.UATBX(IX,IY)
            !!! UBTAX
            New_OcnStat.SST(IX,IY) = &
            New_OcnStat.SST(IX,IY) + (GrdInfo.ReDT) * AllFrc_OcnStat.UBTAX(IX,IY) 
            !!! UATAX
            New_OcnStat.SST(IX,IY) = &
            New_OcnStat.SST(IX,IY) + (GrdInfo.ReDT) * AllFrc_OcnStat.UATAX(IX,IY) 
            
            !!! VATBY
            New_OcnStat.SST(IX,IY) = &
            New_OcnStat.SST(IX,IY) + (GrdInfo.ReDT) * AllFrc_OcnStat.VATBY(IX,IY) 
            !!! VBTAY
            New_OcnStat.SST(IX,IY) = &
            New_OcnStat.SST(IX,IY) + (GrdInfo.ReDT) * AllFrc_OcnStat.VBTAY(IX,IY)
            !!! VATAY
            New_OcnStat.SST(IX,IY) = &
            New_OcnStat.SST(IX,IY) + (GrdInfo.ReDT) * AllFrc_OcnStat.VATAY(IX,IY)
            
            !!! WATBZ
            New_OcnStat.SST(IX,IY) = &
            New_OcnStat.SST(IX,IY) + (GrdInfo.ReDT) * AllFrc_OcnStat.WATBZ(IX,IY)
            !!! WBTAZ
            New_OcnStat.SST(IX,IY) = &
            New_OcnStat.SST(IX,IY) + (GrdInfo.ReDT) * AllFrc_OcnStat.WBTAZ(IX,IY)
            !!! WATAZ
            New_OcnStat.SST(IX,IY) = &
            New_OcnStat.SST(IX,IY) + (GrdInfo.ReDT) * AllFrc_OcnStat.WATAZ(IX,IY)
            
            !!! RESQ
            New_OcnStat.SST(IX,IY) = &
            New_OcnStat.SST(IX,IY) + (GrdInfo.ReDT) * AllFrc_OcnStat.SSTA_RESQ(IX,IY)
            
           

        END DO
     END DO
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !$omp end parallel do
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     
     !!! Sponge Boundary Condition at North/South boundary
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !$omp parallel do schedule (static) &
     !$omp default(shared) &
     !$omp private(IY)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     DO IY = 1, GrdInfo.Sponge_Num
     
        New_OcnStat.SST(1:NX+2,IY) = &
        (1.0d0 - GrdInfo.Sponge_Alp(IY) ) * New_OcnStat.SST(1:NX+2,IY)
        
        New_OcnStat.SST(1:NX+2,NY+2-(IY-1)) = &
        (1.0d0 - GrdInfo.Sponge_Alp(IY) ) * New_OcnStat.SST(1:NX+2,NY+2-(IY-1))
     
     END DO
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !$omp end parallel do
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF ( ( SSTA_Nudging_OnOff == 1 ).AND.(Eq_Direct_Assign == 1) ) THEN
     
     !! Equatorial Nudging (Boundary Condition)
     !$omp parallel do schedule (static) &
     !$omp default(shared) &
     !$omp private(IX,IY)
     DO IX = 1,NX+2
        DO IY = 1,NY+2
     
        IF( (GrdInfo.PLon(IX) >= Nudging_WLon ).AND. ( GrdInfo.PLon(IX) <= Nudging_ELon ) )THEN
     
           IF( (GrdInfo.PLat(IY) >= Nudging_SLat ).AND. ( GrdInfo.PLat(IY) <= Nudging_NLat ) )THEN
        
              New_OcnStat.SST(IX,IY) = New_OcnStat.Nudging_SST(IX,IY) 
           
           END IF
       
        END IF 
     
        END DO
     
     END DO
     !$omp end parallel do
     
     END IF
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     !!! Test nudging pattern 
     !!! New_OcnStat.SST = New_OcnStat.Nudging_SST;
     
     !!! MASK 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !$omp parallel do schedule (static) &
     !$omp default(shared) &
     !$omp private(IX,IY)
     DO IX = 1,NX+2
        DO IY = 1,NY+2
        
        !!! MASK
        New_OcnStat.SST(IX,IY) = New_OcnStat.SST(IX,IY) * GrdInfo.PMask(IX,IY)
        
        END DO
     END DO
     !$omp end parallel do
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     
END SUBROUTINE Add_AllFrc_OceanStat_SSTA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Get_Ekman_Wsrf(GrdInfo,Wsrf,TAUX,TAUY)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Ekman Layer dynamics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! R*(Us) - Beta*Y*(Vs) = (TAUX)/(Rho1*Hm) = Fx
!!! R*(Vs) + Beta*Y*(Us) = (TAUY)/(Rho1*Hm) = Fy
!!!
!!! Us = (R*Fx + Beta*Y*Fy)/(R^2 + (Beta*Y)^2);
!!! Vs = (R*Fy - Beta*Y*Fx)/(R^2 + (Beta*Y)^2);
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TYPE(GridInfo) :: GrdInfo
    real*8,allocatable :: Wsrf(:,:)
    real*8,allocatable :: TAUX(:,:),TAUY(:,:)
    
    real*8,allocatable :: Ekm_Ushr(:,:),Ekm_Vshr(:,:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer*8 :: NX,NY
    integer*8 :: IX,IY
    real*8 :: DX,DY
    
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    
    DX = GrdInfo.ReDX
    DY = GrdInfo.ReDY
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! ekman layer shear
    allocate(Ekm_Ushr(NX+1,NY)); Ekm_Ushr = 0.0d0;
    allocate(Ekm_Vshr(NX,NY+1)); Ekm_Vshr = 0.0d0;
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IX,IY)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO IX = 2,NX 
       DO IY = 1,NY
       
        Ekm_Ushr(IX,IY) = ( Ekm_UVshr_Damp * TAUX(IX,IY) ) + &
        (+0.25d0) *(EarthBeta) * (GrdInfo.UReY(IY)) * &
        ( TAUY(IX-1,IY) + TAUY(IX-1,IY+1) + TAUY(IX,IY) + TAUY(IX,IY+1) )
       
        
        !!! reference mixed layer depth 50.0d0 m
        Ekm_Ushr(IX,IY) = Ekm_Ushr(IX,IY)/( RhoWat * (50.0d0));
        
        
        Ekm_Ushr(IX,IY) = Ekm_Ushr(IX,IY)/ &
        ( Ekm_UVshr_Damp**(2.0d0) + (EarthBeta * GrdInfo.UReY(IY))**(2.0d0) );
        
        
        Ekm_Ushr(IX,IY) = Ekm_Ushr(IX,IY) * (GrdInfo.UMask(IX,IY));
        
       
       END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IX,IY)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO IX = 1,NX
       DO IY = 2,NY
       
       
          Ekm_Vshr(IX,IY) = ( Ekm_UVshr_Damp * TAUY(IX,IY) ) + &
          (-0.25d0) *(EarthBeta) * ( GrdInfo.VReY(IY) ) * &
          ( TAUX(IX,IY-1)  + TAUX(IX,IY) + TAUX(IX+1,IY-1) + TAUX(IX+1,IY) )
          
          
          !!! reference mixed layer depth 50.0d0 m
          Ekm_Vshr(IX,IY) = Ekm_Vshr(IX,IY)/( RhoWat * (50.0d0));
        
        
          Ekm_Vshr(IX,IY) = Ekm_Vshr(IX,IY)/ &
          ( Ekm_UVshr_Damp**(2.0d0) + (EarthBeta * GrdInfo.VReY(IY))**(2.0d0) );
        
        
          Ekm_Vshr(IX,IY) = Ekm_Vshr(IX,IY) * (GrdInfo.VMask(IX,IY));
       
       END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel do schedule (static) &
    !$omp default(shared) &
    !$omp private(IX,IY)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! get Wsrf
    DO IX = 1,NX
       DO IY = 1,NY
       
          Wsrf(IX+1,IY+1) = ( Ekm_Ushr(IX+1,IY) - Ekm_Ushr(IX,IY) )/DX + &
                            ( Ekm_Vshr(IX,IY+1) - Ekm_Vshr(IX,IY) )/DY;
          
          !!! reference mixed layer depth 50.0d0 m
          Wsrf(IX+1,IY+1) = Wsrf(IX+1,IY+1) * (50.0d0); 
    
       END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp end parallel do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    deallocate( Ekm_Ushr, Ekm_Vshr)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE Get_Ekman_Wsrf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ReCover_Surface_DynField (GrdInfo,Clim_OcnStat, TimTick, OcnStat, SrfFlux)
!!! recover real surface physical field 

    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat 
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(SurfaceFlux) :: SrfFlux
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer*8 :: NMode
    integer*8 :: IM 
    
    integer*8 :: NX,NY,NZ,ULNZ
    integer*8 :: IX,IY,IZ
    
    integer*8 :: CMonth !!! the calendar month
    
    real*8 :: DX,DY
    
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    NZ = GrdInfo.NZ
    ULNZ = GrdInfo.ULNZ
    NMode = GrdInfo.NMode
    
    DX = GrdInfo.ReDX
    DY = GrdInfo.ReDY
    
    CMonth = TimTick.CMonth
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! recover SSH 
    OcnStat.SSH = 0.0d0;
    OcnStat.UL_DzDt = 0.0d0;
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO IM = 1,NMode
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       DO IX = 1,NX
          DO IY = 1,NY
          
             !!! be careful EigMode(IX,1,IM,CMonth)
             OcnStat.SSH(IX+1,IY+1) = OcnStat.SSH(IX+1,IY+1) + &
             ( Clim_OcnStat.EigMode(IX+1,1,IM) / EarthG ) * OcnStat.Co_Pres(IX+1,IY+1,IM); 
    
          END DO
       END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    END DO
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! recover Usrf 
    OcnStat.Usrf = 0.0d0;
    OcnStat.UL_Ucur = 0.0d0;
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO IM = 1,NMode
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       DO IX = 1,NX+1
          DO IY = 1,NY
          
              OcnStat.Usrf(IX,IY) = OcnStat.Usrf(IX,IY) + &
             ( Clim_OcnStat.EigMode(IX+1,1,IM) ) * OcnStat.Co_Ucur(IX,IY,IM); 
          
          
          END DO
       END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Recover Vsrf 
    OcnStat.Vsrf = 0.0d0;
    OcnStat.UL_Vcur = 0.0d0;
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO IM = 1,NMode
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       DO IX = 1,NX
          DO IY = 1,NY+1
          
              OcnStat.Vsrf(IX,IY) = OcnStat.Vsrf(IX,IY) + &
             ( Clim_OcnStat.EigMode(IX+1,1,IM) ) * OcnStat.Co_Vcur(IX,IY,IM); 
          
          END DO
       END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    END DO
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Ekman dynamics
    call Get_Ekman_Wsrf(GrdInfo, OcnStat.Wsrf, SrfFlux.TAUX, SrfFlux.TAUY)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!! OcnStat.SSH = SrfFlux.SSH;    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! OceanStat.TSUB dynamics
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (Solve_TSUB_AI_OnOff == 1 ) THEN
    
        CALL Get_TSUB_AI(GrdInfo, TimTick, OcnStat, Clim_OcnStat)
        
        !!! Fotran-Torch-Adapter BUG (need activated)
        IF( TimTick.Run_AbsIT == (TimTick.Res_AbsIT + 1) )THEN
            CALL Get_TSUB_AI(GrdInfo, TimTick, OcnStat, Clim_OcnStat)
        END IF 
        
    ELSE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !DO IX = 1,NX
    !   DO IY = 1,NY
    !   
    !      !!!! control thermocline feedback 
    !      
    !      OcnStat.TSUB(IX+1,IY+1) = Get_TSUB ( SSHA = OcnStat.SSH(IX+1,IY+1), &
    !                                           SSTA = OcnStat.SST(IX+1,IY+1), &
    !                                           TCDB = Clim_OcnStat.TCDBar(IX+1,IY+1,CMonth) );
    !   
    !   END DO
    !END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
END SUBROUTINE ReCover_Surface_DynField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ReCover_UpperLayer_DynField (GrdInfo,Clim_OcnStat, TimTick, OcnStat)
!!! recover real UpperLayer physical field (optional)

    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat 
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer*8 :: NMode
    integer*8 :: IM 
    
    integer*8 :: NX,NY,NZ,ULNZ
    integer*8 :: IX,IY,IZ
    
    integer*8 :: CMonth !!! the calendar month
    
    real*8 :: DX,DY
    
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    NZ = GrdInfo.NZ
    ULNZ = GrdInfo.ULNZ
    NMode = GrdInfo.NMode
    
    DX = GrdInfo.ReDX
    DY = GrdInfo.ReDY
    
    CMonth = TimTick.CMonth
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! recover SSH & UL_DzDt 
    OcnStat.SSH = 0.0d0;
    OcnStat.UL_DzDt = 0.0d0;
    
    DO IM = 1,NMode
       DO IX = 1,NX
          DO IY = 1,NY
          
             !!! using IZ_EigMode
             DO IZ = 1,ULNZ
             
                OcnStat.UL_DzDt(IX+1,IY+1,IZ) = OcnStat.UL_DzDt(IX+1,IY+1,IZ) + &
                ( Clim_OcnStat.IZ_EigMode(IX+1,IZ,IM) ) * &
                ( (1.0d0/DX) * ( OcnStat.Co_Ucur(IX+1,IY,IM) - OcnStat.Co_Ucur(IX,IY,IM) ) + &
                  (1.0d0/DY) * ( OcnStat.Co_Vcur(IX,IY+1,IM) - OcnStat.Co_Vcur(IX,IY,IM) ) );
                
             END DO
             
          END DO
       END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Recover Usrf & Ucur & 
    OcnStat.Usrf = 0.0d0;
    OcnStat.UL_Ucur = 0.0d0;
    
    DO IM = 1,NMode
       DO IX = 1,NX+1
          DO IY = 1,NY
             DO IZ = 1,ULNZ
             
                OcnStat.UL_Ucur(IX,IY,IZ) = OcnStat.UL_Ucur(IX,IY,IZ) + &
                ( Clim_OcnStat.EigMode(IX+1,IZ,IM) ) * OcnStat.Co_Ucur(IX,IY,IM); 
             
             END DO
          END DO
       END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Recover Vsrf & Vcur & 
    OcnStat.Vsrf = 0.0d0;
    OcnStat.UL_Vcur = 0.0d0;
    
    DO IM = 1,NMode
       DO IX = 1,NX
          DO IY = 1,NY+1
          
             DO IZ = 1,ULNZ
             
                OcnStat.UL_Vcur(IX,IY,IZ) = OcnStat.UL_Vcur(IX,IY,IZ) + &
                ( Clim_OcnStat.EigMode(IX+1,IZ,IM) ) * OcnStat.Co_Vcur(IX,IY,IM); 
             
             END DO
          
          END DO
       END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE ReCover_UpperLayer_DynField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Get_NinoIndex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Get_NinoIndex (GrdInfo,OcnStat, SrfFlux)
    
    TYPE(GridInfo) :: GrdInfo 
    TYPE(OceanStat) :: OcnStat
    TYPE(SurfaceFlux) :: SrfFlux
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel sections
    !$omp section
    
    OcnStat.SSHA_N4 = sum(OcnStat.SSH( GrdInfo.Nino4_Bd.WInd : GrdInfo.Nino4_Bd.EInd, &
                                       GrdInfo.Nino4_Bd.SInd : GrdInfo.Nino4_Bd.NInd ) );
    OcnStat.SSHA_N4 = OcnStat.SSHA_N4/dble(GrdInfo.Nino4_PointNum);
    
     
    OcnStat.SSHA_N34 = sum(OcnStat.SSH( GrdInfo.Nino3p4_Bd.WInd : GrdInfo.Nino3p4_Bd.EInd, &
                                       GrdInfo.Nino3p4_Bd.SInd : GrdInfo.Nino3p4_Bd.NInd ) );
    OcnStat.SSHA_N34 = OcnStat.SSHA_N34/dble(GrdInfo.Nino3p4_PointNum);
    
    
    OcnStat.SSHA_N3 = sum(OcnStat.SSH( GrdInfo.Nino3_Bd.WInd : GrdInfo.Nino3_Bd.EInd, &
                                       GrdInfo.Nino3_Bd.SInd : GrdInfo.Nino3_Bd.NInd ) );
    OcnStat.SSHA_N3 = OcnStat.SSHA_N3/dble(GrdInfo.Nino3_PointNum);
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !$omp section
    
    OcnStat.SSTA_N4 = sum(OcnStat.SST( GrdInfo.Nino4_Bd.WInd : GrdInfo.Nino4_Bd.EInd, &
                                       GrdInfo.Nino4_Bd.SInd : GrdInfo.Nino4_Bd.NInd ) );
    OcnStat.SSTA_N4 = OcnStat.SSTA_N4/dble(GrdInfo.Nino4_PointNum);
    
     
    OcnStat.SSTA_N34 = sum(OcnStat.SST( GrdInfo.Nino3p4_Bd.WInd : GrdInfo.Nino3p4_Bd.EInd, &
                                       GrdInfo.Nino3p4_Bd.SInd : GrdInfo.Nino3p4_Bd.NInd ) );
    OcnStat.SSTA_N34 = OcnStat.SSTA_N34/dble(GrdInfo.Nino3p4_PointNum);
    
    
    OcnStat.SSTA_N3 = sum(OcnStat.SST( GrdInfo.Nino3_Bd.WInd : GrdInfo.Nino3_Bd.EInd, &
                                       GrdInfo.Nino3_Bd.SInd : GrdInfo.Nino3_Bd.NInd ) );
    OcnStat.SSTA_N3 = OcnStat.SSTA_N3/dble(GrdInfo.Nino3_PointNum);
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !$omp section
    
    OcnStat.Eq_SSTA_N4 = sum(OcnStat.SST( GrdInfo.Eq_Nino4_Bd.WInd : GrdInfo.Eq_Nino4_Bd.EInd, &
                                         GrdInfo.Eq_Nino4_Bd.SInd : GrdInfo.Eq_Nino4_Bd.NInd ) );
    OcnStat.Eq_SSTA_N4 = OcnStat.Eq_SSTA_N4/dble(GrdInfo.Eq_Nino4_PointNum);
    
     
    OcnStat.Eq_SSTA_N34 = sum(OcnStat.SST( GrdInfo.Eq_Nino3p4_Bd.WInd : GrdInfo.Eq_Nino3p4_Bd.EInd, &
                                              GrdInfo.Eq_Nino3p4_Bd.SInd : GrdInfo.Eq_Nino3p4_Bd.NInd ) );
    OcnStat.Eq_SSTA_N34 = OcnStat.Eq_SSTA_N34/dble(GrdInfo.Eq_Nino3p4_PointNum);
    
    
    OcnStat.Eq_SSTA_N3 = sum(OcnStat.SST( GrdInfo.Eq_Nino3_Bd.WInd : GrdInfo.Eq_Nino3_Bd.EInd, &
                                            GrdInfo.Eq_Nino3_Bd.SInd : GrdInfo.Eq_Nino3_Bd.NInd ) );
    OcnStat.Eq_SSTA_N3 = OcnStat.Eq_SSTA_N3/dble(GrdInfo.Eq_Nino3_PointNum);
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !$omp section
    OcnStat.TSUB_N4 = sum(OcnStat.TSUB( GrdInfo.Nino4_Bd.WInd : GrdInfo.Nino4_Bd.EInd, &
                                       GrdInfo.Nino4_Bd.SInd : GrdInfo.Nino4_Bd.NInd ) );
    OcnStat.TSUB_N4 = OcnStat.TSUB_N4/dble(GrdInfo.Nino4_PointNum);
    
     
    OcnStat.TSUB_N34 = sum(OcnStat.TSUB( GrdInfo.Nino3p4_Bd.WInd : GrdInfo.Nino3p4_Bd.EInd, &
                                       GrdInfo.Nino3p4_Bd.SInd : GrdInfo.Nino3p4_Bd.NInd ) );
    OcnStat.TSUB_N34 = OcnStat.TSUB_N34/dble(GrdInfo.Nino3p4_PointNum);
    
    
    OcnStat.TSUB_N3 = sum(OcnStat.TSUB( GrdInfo.Nino3_Bd.WInd : GrdInfo.Nino3_Bd.EInd, &
                                       GrdInfo.Nino3_Bd.SInd : GrdInfo.Nino3_Bd.NInd ) );
    OcnStat.TSUB_N3 = OcnStat.TSUB_N3/dble(GrdInfo.Nino3_PointNum);
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! be careful for Arakawa-C Grid
    !!!          U(IX,IY)
    !!! V(IX,IY) P(IX+1,IY+1) V(IX,IY+1)
    !!!          U(IX+1,IY)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !$omp section
    
    OcnStat.UA_N4 = sum(OcnStat.Usrf( GrdInfo.Nino4_Bd.WInd-1 : GrdInfo.Nino4_Bd.EInd-1, &
                                      GrdInfo.Nino4_Bd.SInd-1 : GrdInfo.Nino4_Bd.NInd-1 ) );
    OcnStat.UA_N4 = OcnStat.UA_N4/dble(GrdInfo.Nino4_PointNum);
    
    OcnStat.UA_N34 = sum(OcnStat.Usrf( GrdInfo.Nino3p4_Bd.WInd-1 : GrdInfo.Nino3p4_Bd.EInd-1, &
                                       GrdInfo.Nino3p4_Bd.SInd-1 : GrdInfo.Nino3p4_Bd.NInd-1 ) );
    OcnStat.UA_N34 = OcnStat.UA_N34/dble(GrdInfo.Nino3p4_PointNum);
    
    OcnStat.UA_N3 = sum(OcnStat.Usrf( GrdInfo.Nino3_Bd.WInd-1 : GrdInfo.Nino3_Bd.EInd-1, &
                                      GrdInfo.Nino3_Bd.SInd-1 : GrdInfo.Nino3_Bd.NInd-1 ) );
    OcnStat.UA_N3 = OcnStat.UA_N3/dble(GrdInfo.Nino3_PointNum);
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !$omp section
    
    !!! TAUX N4,N34,N3 series
    
    SrfFlux.TAUX_N4 = sum(SrfFlux.TAUX( GrdInfo.Nino4_Bd.WInd-1 : GrdInfo.Nino4_Bd.EInd-1, &
                                        GrdInfo.Nino4_Bd.SInd-1 : GrdInfo.Nino4_Bd.NInd-1 ) );
    SrfFlux.TAUX_N4 = SrfFlux.TAUX_N4/dble(GrdInfo.Nino4_PointNum);
    
    SrfFlux.TAUX_N34 = sum(SrfFlux.TAUX( GrdInfo.Nino3p4_Bd.WInd-1 : GrdInfo.Nino3p4_Bd.EInd-1, &
                                         GrdInfo.Nino3p4_Bd.SInd-1 : GrdInfo.Nino3p4_Bd.NInd-1 ) );
    SrfFlux.TAUX_N34 = SrfFlux.TAUX_N34/dble(GrdInfo.Nino3p4_PointNum);
    
    SrfFlux.TAUX_N3 = sum(SrfFlux.TAUX( GrdInfo.Nino3_Bd.WInd-1 : GrdInfo.Nino3_Bd.EInd-1, &
                                        GrdInfo.Nino3_Bd.SInd-1 : GrdInfo.Nino3_Bd.NInd-1 ) );
    SrfFlux.TAUX_N3 = SrfFlux.TAUX_N3/dble(GrdInfo.Nino3_PointNum);
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !$omp section
    
    !!! TAUX N4,N34,N3 test series
    
    SrfFlux.TAUX_N4_Test = sum(SrfFlux.TAUX_Test( GrdInfo.Nino4_Bd.WInd-1 : GrdInfo.Nino4_Bd.EInd-1, &
                                                  GrdInfo.Nino4_Bd.SInd-1 : GrdInfo.Nino4_Bd.NInd-1 ) );
    SrfFlux.TAUX_N4_Test = SrfFlux.TAUX_N4_Test/dble(GrdInfo.Nino4_PointNum);
    
    SrfFlux.TAUX_N34_Test = sum(SrfFlux.TAUX_Test( GrdInfo.Nino3p4_Bd.WInd-1 : GrdInfo.Nino3p4_Bd.EInd-1, &
                                                   GrdInfo.Nino3p4_Bd.SInd-1 : GrdInfo.Nino3p4_Bd.NInd-1 ) );
    SrfFlux.TAUX_N34_Test = SrfFlux.TAUX_N34_Test/dble(GrdInfo.Nino3p4_PointNum);
    
    SrfFlux.TAUX_N3_Test = sum(SrfFlux.TAUX_Test( GrdInfo.Nino3_Bd.WInd-1 : GrdInfo.Nino3_Bd.EInd-1, &
                                                  GrdInfo.Nino3_Bd.SInd-1 : GrdInfo.Nino3_Bd.NInd-1 ) );
    SrfFlux.TAUX_N3_Test = SrfFlux.TAUX_N3_Test/dble(GrdInfo.Nino3_PointNum);
    
    !$omp end parallel sections 
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
END SUBROUTINE Get_NinoIndex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Get_MonMe_Stat(TimTick, GrdInfo, OcnStat, Clim_OcnStat, AllFrc_OcnStat,SrfFlux)
!!! get monthly mean stat

    TYPE(TimeTick) :: TimTick
    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat 
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(AllFrc_OceanStat) :: AllFrc_OcnStat
    TYPE(SurfaceFlux) :: SrfFlux
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer*8 :: NMode
    integer*8 :: IM 
    
    integer*8 :: NX,NY,NZ,ULNZ
    integer*8 :: IX,IY,IZ
    
    integer*8 :: Arrive_CMonth_FLAG !!! Arrive_FLAG: arrive at the start of month
    integer*8 :: Reach_CMonth_FLAG  !!! Reach_FLAG:  arrive at the end of month
    
    integer*8 :: CYear, CMonth

    
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    NZ = GrdInfo.NZ
    NMode = GrdInfo.NMode

    Arrive_CMonth_FLAG = TimTick.Arrive_CMonth_FLAG
    Reach_CMonth_FLAG  = TimTick.Reach_CMonth_FLAG
    
    CYear  = TimTick.CYear
    CMonth = TimTick.CMonth
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !$omp parallel sections
    
    !!! MonMe OcnStat
    !$omp section
    call Update_MonMe_Var (OcnStat.Usrf, OcnStat.MonMe_Usrf, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    !$omp section
    call Update_MonMe_Var (OcnStat.Vsrf, OcnStat.MonMe_Vsrf, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    !$omp section
    call Update_MonMe_Var (OcnStat.Wsrf, OcnStat.MonMe_Wsrf, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    !$omp section
    call Update_MonMe_Var (OcnStat.SSH, OcnStat.MonMe_SSH, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    !$omp section
    call Update_MonMe_Var (OcnStat.SST, OcnStat.MonMe_SST, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    !$omp section
    call Update_MonMe_Var (OcnStat.TSUB, OcnStat.MonMe_TSUB, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    call Update_MonMe_Var (OcnStat.SSTA_N4, OcnStat.MonMe_SSTA_N4, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    call Update_MonMe_Var (OcnStat.SSTA_N34, OcnStat.MonMe_SSTA_N34, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    call Update_MonMe_Var (OcnStat.SSTA_N3, OcnStat.MonMe_SSTA_N3, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    
    !!! MonMe AllFrc_OcnStat
    !$omp section
    call Update_MonMe_Var (AllFrc_OcnStat.ZAFK, AllFrc_OcnStat.MonMe_ZAFK, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    !$omp section
    call Update_MonMe_Var (AllFrc_OcnStat.THFK, AllFrc_OcnStat.MonMe_THFK, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    !$omp section
    call Update_MonMe_Var (AllFrc_OcnStat.RESFK, AllFrc_OcnStat.MonMe_RESFK, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    !$omp section
    call Update_MonMe_Var (AllFrc_OcnStat.NUDGE, AllFrc_OcnStat.MonMe_NUDGE, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    
    !!! MonMe SrfFlux
    !$omp section
    call Update_MonMe_Var (SrfFlux.TAUX, SrfFlux.MonMe_TAUX, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    !$omp section
    call Update_MonMe_Var (SrfFlux.TAUY, SrfFlux.MonMe_TAUY, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    !$omp section
    call Update_MonMe_Var (SrfFlux.TAUX_Test, SrfFlux.MonMe_TAUX_Test, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    !$omp section
    call Update_MonMe_Var (SrfFlux.TAUY_Test, SrfFlux.MonMe_TAUY_Test, &
         Arrive_CMonth_FLAG, Reach_CMonth_FLAG, Leap_Calendar_OnOff, CYear, CMonth, Day_NT)
    
    !$omp end parallel sections 
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE Get_MonMe_Stat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OverFlow_Judge (GrdInfo,OcnStat,SrfFlux)

    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat
    TYPE(SurfaceFlux) :: SrfFlux
    
    integer*8 :: NX,NY,NZ
    integer*8 :: IX,IY,IZ,IM
    
    NX = GrdInfo.NX
    NY = GrdInfo.NY
    
    IF ( abs(OcnStat.SSHA_N3) >= 1000.0d0 ) THEN 
        write(*,*),'Numerical Value Overflow, Program Stopped'
        STOP
    END IF
    
    IF ( abs(OcnStat.SSTA_N3) >= 1000.0d0 ) THEN 
        write(*,*),'Numerical Value Overflow, Program Stopped'
        STOP
    END IF

END SUBROUTINE OverFlow_Judge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OneStep_Forward &
           (GrdInfo,Clim_OcnStat, TimTick, OcnStat, Tmp_OcnStat, AllFrc_OcnStat, SrfFlux)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TYPE(GridInfo) :: GrdInfo
    TYPE(OceanStat) :: OcnStat, Tmp_OcnStat 
    TYPE(Clim_OceanStat) :: Clim_OcnStat
    TYPE(TimeTick) :: TimTick
    TYPE(AllFrc_OceanStat) :: AllFrc_OcnStat
    TYPE(SurfaceFlux) :: SrfFlux
    
    integer*8 :: IX,IY,IT
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! Mastuno two-step Scheme for dynamic field
    
    !!!! F(Un)
    call Get_AllFrc_OceanStat_MulMode &
         (GrdInfo,Clim_OcnStat, TimTick, OcnStat, AllFrc_OcnStat,SrfFlux)
    
    !!! U*  = Un + DT * F(Un)
    call Add_AllFrc_OceanStat_MulMode &
         (GrdInfo,Clim_OcnStat, TimTick, OcnStat, Tmp_OcnStat, AllFrc_OcnStat)
    
    !!!! F(U*) 
    call Get_AllFrc_OceanStat_MulMode &
         (GrdInfo,Clim_OcnStat, TimTick, Tmp_OcnStat, AllFrc_OcnStat,SrfFlux)
    
    !!! Un+1 = Un + DT * F(U*)
    call Add_AllFrc_OceanStat_MulMode &
         (GrdInfo,Clim_OcnStat, TimTick, OcnStat, OcnStat, AllFrc_OcnStat)
    
    !!! Recover Surface Dynamic Field (Usrf, Vsrf, Wsrf)
    call ReCover_Surface_DynField (GrdInfo,Clim_OcnStat, TimTick, OcnStat, SrfFlux)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    IF (Solve_SSTA_OnOff == 1 ) THEN
    
    !!!! Euler forward for SSTA tendency equation
    !!! F(Un)
    call Get_AllFrc_OceanStat_SSTA &
         (GrdInfo,Clim_OcnStat, TimTick, OcnStat, AllFrc_OcnStat,SrfFlux)
    
    !!! Un+1 = Un + DT * F(U)
    call Add_AllFrc_OceanStat_SSTA &
         (GrdInfo,Clim_OcnStat, TimTick, OcnStat, OcnStat, AllFrc_OcnStat)
    
    END IF
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!  Nino Index
    call Get_NinoIndex (GrdInfo,OcnStat,SrfFlux)
   
    !!! OverFlow Judge
    call OverFlow_Judge (GrdInfo,OcnStat,SrfFlux)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!! Recover UpperLayer Physical Field in (X,Y,Z,T) space
    !!! Only Recover at OutPut Time
    IF (OutPut_UpperLayer_OnOff == 1) THEN
    
    IF ( (TimTick.Run_AbsIT>=1) .AND. (mod(TimTick.Run_AbsIT, (Day_NT)*(PutData_N_Day)) == 0) )THEN
    
        call ReCover_UpperLayer_DynField (GrdInfo,Clim_OcnStat, TimTick, OcnStat)
        
    END IF
        
    END IF !!! OutPut_UpperLayer_OnOff
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE OneStep_Forward
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE MO_OCEANSOLVER
