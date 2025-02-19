MODULE MO_NCTOOLS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NetCDF ToolBox
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    USE TypeSizes
    USE NetCDF
    
    IMPLICIT NONE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! OVERLOAD
    INTERFACE NC_Read_Mul_IT
        
        module procedure NC_Read_1Dto1D
        module procedure NC_Read_2Dto2D
        module procedure NC_Read_3Dto3D
        module procedure NC_Read_4Dto4D
        module procedure NC_Read_5Dto5D
        
    END INTERFACE NC_Read_Mul_IT
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTERFACE NC_Read_One_IT
    
        module procedure NC_Read_1Dto0D
        module procedure NC_Read_2Dto1D
        module procedure NC_Read_3Dto2D
        module procedure NC_Read_4Dto3D
        module procedure NC_Read_5Dto4D
    
    END INTERFACE NC_Read_One_IT
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTERFACE NC_Create_File
       
        module procedure NC_Create_1D_File
        module procedure NC_Create_2D_File
        module procedure NC_Create_3D_File
        module procedure NC_Create_4D_File
        module procedure NC_Create_5D_File 
    
    END INTERFACE NC_Create_File
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTERFACE NC_Def_Var
       
        module procedure NC_Def_1D_Var
        module procedure NC_Def_2D_Var
        module procedure NC_Def_3D_Var
        module procedure NC_Def_4D_Var
        module procedure NC_Def_5D_Var
    
    END INTERFACE NC_Def_Var
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTERFACE NC_Write_Mul_IT
    
        module procedure NC_Write_1Dto1D
        module procedure NC_Write_2Dto2D
        module procedure NC_Write_3Dto3D
        module procedure NC_Write_4Dto4D
        module procedure NC_Write_5Dto5D
    
    END INTERFACE NC_Write_Mul_IT
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTERFACE NC_Write_One_IT
    
        module procedure NC_Write_0Dto1D
        module procedure NC_Write_1Dto2D
        module procedure NC_Write_2Dto3D
        module procedure NC_Write_3Dto4D
        module procedure NC_Write_4Dto5D
    
    END INTERFACE NC_Write_One_IT
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    CONTAINS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_check
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Check(status)
    
        integer*4 :: status
    
        if(status /= nf90_noerr) then 
          print *, trim(nf90_strerror(status))
        end if
        
    END SUBROUTINE NC_Check  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!
    !!! NC_Read_Mul_IT
    !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Read_1Dto1D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Read_1Dto1D(FPath,FName,VarName,PutArray,StartIT,CountT,FLAG)
    
        character(len=*) :: FPath,FName,VarName
        real*8,allocatable :: PutArray(:)
        
        integer*4 :: ncID,varID
        integer*4, dimension(nf90_max_var_dims) :: dimIDs
        integer*4 :: Dim01
        integer*4 :: StartIT,CountT
        character :: FLAG !!! FLAG == 'A':allocate, FLAG == 'N',noallocate
        
        call NC_Check(nf90_open(trim(FPath)//trim(FName),nf90_NoWrite,ncID));
        
        call NC_Check(nf90_inq_varid(ncID,trim(VarName),varID));
        
        call NC_Check(nf90_Inquire_Variable(ncID,varID,dimids = dimIDs))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(1), len=Dim01))
        
         IF( StartIT<= 0 )THEN
            IF( FLAG == 'A')THEN
            allocate(PutArray(Dim01));
            END IF
            call NC_Check(nf90_get_var(ncID, varID, PutArray));
        ELSE
            IF( FLAG == 'A')THEN
            allocate(PutArray(CountT));
            END IF    
            call NC_Check(nf90_get_var ( ncID, varID,PutArray,(/int4(StartIT)/),(/int4(CountT)/) ) );
        END IF
        
        call NC_Check(nf90_close(ncID))
        
        write(*,*) 'Reading '//trim(VarName)//' from '//trim(FPath)//trim(FName)
        write(*,*),''
        
    END SUBROUTINE NC_Read_1Dto1D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Read_2Dto2D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Read_2Dto2D(FPath,FName,VarName,PutArray,StartIT,CountT,FLAG)
    
        character(len=*) :: FPath,FName,VarName
        real*8,allocatable :: PutArray(:,:)
        
        integer*4 :: ncID,varID
        integer*4, dimension(nf90_max_var_dims) :: dimIDs
        integer*4 :: Dim01,Dim02
        integer*4 :: StartIT, CountT
        character :: FLAG !!! FLAG == 'A':allocate, FLAG == 'N',noallocate
        
        call NC_Check(nf90_open(trim(FPath)//trim(FName),nf90_NoWrite,ncID));
        
        call NC_Check(nf90_inq_varid(ncID,trim(VarName),varID));
        
        call NC_Check(nf90_Inquire_Variable(ncID,varID,dimids = dimIDs))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(1), len=Dim01))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(2), len=Dim02))
        
        IF( StartIT<= 0 )THEN
            IF(FLAG == 'A')THEN
            allocate(PutArray(Dim01,Dim02));
            END IF
            
            call NC_Check(nf90_get_var(ncID, varID, PutArray));
        ELSE
            IF( FLAG == 'A')THEN
            allocate(PutArray(Dim01,CountT));
            END IF
            
            call NC_Check(nf90_get_var ( ncID, varID,PutArray,&
                          (/1,int4(StartIT)/),(/Dim01,int4(CountT)/) ) );
        END IF
        
        call NC_Check(nf90_close(ncID))
        
        write(*,*) 'Reading '//trim(VarName)//' from '//trim(FPath)//trim(FName)
        write(*,*),''
        
    END SUBROUTINE NC_Read_2Dto2D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Read_3Dto3D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Read_3Dto3D(FPath,FName,VarName,PutArray,StartIT,CountT,FLAG)
    
        character(len=*) :: FPath,FName,VarName
        real*8,allocatable :: PutArray(:,:,:)
        
        integer*4 :: ncID,varID
        integer*4, dimension(nf90_max_var_dims) :: dimIDs
        integer*4 :: Dim01,Dim02,Dim03
        integer*4 :: StartIT, CountT
        character :: FLAG !!! FLAG == 'A':allocate, FLAG == 'N',noallocate
        
        call NC_Check(nf90_open(trim(FPath)//trim(FName),nf90_NoWrite,ncID));
        
        call NC_Check(nf90_inq_varid(ncID,trim(VarName),varID));
        
        call NC_Check(nf90_Inquire_Variable(ncID,varID,dimids = dimIDs))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(1), len=Dim01))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(2), len=Dim02))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(3), len=Dim03))
        
        IF( StartIT<= 0 )THEN
            IF( FLAG == 'A')THEN
            allocate(PutArray(Dim01,Dim02,Dim03));
            END IF
            
            call NC_Check(nf90_get_var(ncID, varID, PutArray));
        ELSE
            IF( FLAG == 'A')THEN
            allocate(PutArray(Dim01,Dim02,CountT));
            END IF
           
            call NC_Check(nf90_get_var ( ncID, varID,PutArray,&
                          (/1,1,int4(StartIT)/),(/Dim01,Dim02,int4(CountT)/) ) );
        END IF
        
        call NC_Check(nf90_close(ncID))
        
        write(*,*) 'Reading '//trim(VarName)//' from '//trim(FPath)//trim(FName)
        write(*,*),''
        
    END SUBROUTINE NC_Read_3Dto3D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Read_4Dto4D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Read_4Dto4D(FPath,FName,VarName,PutArray,StartIT,CountT,FLAG)
    
        character(len=*) :: FPath,FName,VarName
        real*8,allocatable :: PutArray(:,:,:,:)
        
        integer*4 :: ncID,varID
        integer*4, dimension(nf90_max_var_dims) :: dimIDs
        integer*4 :: Dim01,Dim02,Dim03,Dim04
        integer*4 :: StartIT, CountT
        character :: FLAG !!! FLAG == 'A':allocate, FLAG == 'N',noallocate
        
        call NC_Check(nf90_open(trim(FPath)//trim(FName),nf90_NoWrite,ncID));
        
        call NC_Check(nf90_inq_varid(ncID,trim(VarName),varID));
        
        call NC_Check(nf90_Inquire_Variable(ncID,varID,dimids = dimIDs))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(1), len=Dim01))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(2), len=Dim02))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(3), len=Dim03))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(4), len=Dim04))
        
       
        IF( StartIT<= 0 )THEN
            IF( FLAG == 'A')THEN
            allocate(PutArray(Dim01,Dim02,Dim03,Dim04));
            END IF
            call NC_Check(nf90_get_var(ncID, varID, PutArray));
        ELSE
            IF( FLAG == 'A')THEN
            allocate(PutArray(Dim01,Dim02,Dim03,CountT));
            END IF
            call NC_Check(nf90_get_var ( ncID, varID,PutArray,&
                          (/1,1,1,int4(StartIT)/),(/Dim01,Dim02,Dim03,int4(CountT)/) ) );
        END IF
        
        call NC_Check(nf90_close(ncID))
        
        write(*,*) 'Reading '//trim(VarName)//' from '//trim(FPath)//trim(FName)
        write(*,*),''
        
    END SUBROUTINE NC_Read_4Dto4D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Read_4Dto4D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Read_5Dto5D(FPath,FName,VarName,PutArray,StartIT,CountT,FLAG)
    
        character(len=*) :: FPath,FName,VarName
        real*8,allocatable :: PutArray(:,:,:,:,:)
        
        integer*4 :: ncID,varID
        integer*4, dimension(nf90_max_var_dims) :: dimIDs
        integer*4 :: Dim01,Dim02,Dim03,Dim04,Dim05
        integer*4 :: StartIT, CountT
        character :: FLAG !!! FLAG == 'A':allocate, FLAG == 'N',noallocate
        
        call NC_Check(nf90_open(trim(FPath)//trim(FName),nf90_NoWrite,ncID));
        
        call NC_Check(nf90_inq_varid(ncID,trim(VarName),varID));
        
        call NC_Check(nf90_Inquire_Variable(ncID,varID,dimids = dimIDs))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(1), len=Dim01))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(2), len=Dim02))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(3), len=Dim03))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(4), len=Dim04))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(5), len=Dim05))
        
       
        IF( StartIT<= 0 )THEN
            IF( FLAG == 'A')THEN
            allocate(PutArray(Dim01,Dim02,Dim03,Dim04,Dim05));
            END IF
            call NC_Check(nf90_get_var(ncID, varID, PutArray));
        ELSE
            IF( FLAG == 'A')THEN
            allocate(PutArray(Dim01,Dim02,Dim03,Dim04,CountT));
            END IF
            call NC_Check(nf90_get_var ( ncID, varID,PutArray,&
                          (/1,1,1,1,int4(StartIT)/),(/Dim01,Dim02,Dim03,Dim04,int4(CountT)/) ) );
        END IF
        
        call NC_Check(nf90_close(ncID))
        
        write(*,*) 'Reading '//trim(VarName)//' from '//trim(FPath)//trim(FName)
        write(*,*),''
        
    END SUBROUTINE NC_Read_5Dto5D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!
    !!! NC_Read_One_IT
    !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Read_1Dto0D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Read_1Dto0D(FPath,FName,VarName,PutArray,TimeInd,FLAG)
    
        character(len=*) :: FPath,FName,VarName
        real*8 :: PutArray
        
        integer*4 :: ncID,varID
        integer*4, dimension(nf90_max_var_dims) :: dimIDs
        integer*4 :: Dim01
        integer*4 :: TimeInd
        character :: FLAG !!! FLAG == 'A':allocate, FLAG == 'N',noallocate
        
        call NC_Check(nf90_open(trim(FPath)//trim(FName),nf90_NoWrite,ncID));
        
        call NC_Check(nf90_inq_varid(ncID,trim(VarName),varID));
        
        call NC_Check(nf90_Inquire_Variable(ncID,varID,dimids = dimIDs))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(1), len=Dim01))
        
        call NC_Check( nf90_get_var(ncID, varID, PutArray,(/int4(TimeInd)/)) );
        
        call NC_Check(nf90_close(ncID))
        
        write(*,*) 'Reading '//trim(VarName)//' from '//trim(FPath)//trim(FName)
        write(*,*),''
        
    END SUBROUTINE NC_Read_1Dto0D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Read_2Dto1D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Read_2Dto1D(FPath,FName,VarName,PutArray,TimeInd,FLAG)
    
        character(len=*) :: FPath,FName,VarName
        real*8,allocatable :: PutArray(:)
        
        integer*4 :: ncID,varID
        integer*4, dimension(nf90_max_var_dims) :: dimIDs
        integer*4 :: Dim01,Dim02
        integer*4 :: TimeInd
        character :: FLAG !!! FLAG == 'A':allocate, FLAG == 'N',noallocate
        
        call NC_Check(nf90_open(trim(FPath)//trim(FName),nf90_NoWrite,ncID));
        
        call NC_Check(nf90_inq_varid(ncID,trim(VarName),varID));
        
        call NC_Check(nf90_Inquire_Variable(ncID,varID,dimids = dimIDs))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(1), len=Dim01))
        
        IF( FLAG == 'A')THEN 
        allocate(PutArray(Dim01));
        END IF
        
        call NC_Check( nf90_get_var(ncID, varID, PutArray,(/1,int4(TimeInd)/)) );
        
        call NC_Check(nf90_close(ncID))
        
        write(*,*) 'Reading '//trim(VarName)//' from '//trim(FPath)//trim(FName)
        write(*,*),''
        
    END SUBROUTINE NC_Read_2Dto1D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Read_3Dto2D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Read_3Dto2D(FPath,FName,VarName,PutArray,TimeInd,FLAG)
    
        character(len=*) :: FPath,FName,VarName
        real*8,allocatable :: PutArray(:,:)
        
        integer*4 :: ncID,varID
        integer*4, dimension(nf90_max_var_dims) :: dimIDs
        integer*4 :: Dim01,Dim02
        integer*4 :: TimeInd
        character :: FLAG !!! FLAG == 'A':allocate, FLAG == 'N',noallocate
        
        call NC_Check(nf90_open(trim(FPath)//trim(FName),nf90_NoWrite,ncID));
        
        call NC_Check(nf90_inq_varid(ncID,trim(VarName),varID));
        
        call NC_Check(nf90_Inquire_Variable(ncID,varID,dimids = dimIDs))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(1), len=Dim01))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(2), len=Dim02))
        
        
        IF( FLAG == 'A')THEN
        allocate(PutArray(Dim01,Dim02));
        END IF
        
        call NC_Check(nf90_get_var(ncID, varID, PutArray,(/1,1,int4(TimeInd)/)));
        
        call NC_Check(nf90_close(ncID))
        
        write(*,*) 'Reading '//trim(VarName)//' from '//trim(FPath)//trim(FName)
        write(*,*),''
        
    END SUBROUTINE NC_Read_3Dto2D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Read_4Dto3D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Read_4Dto3D(FPath,FName,VarName,PutArray,TimeInd,FLAG)
    
        character(len=*) :: FPath,FName,VarName
        real*8,allocatable :: PutArray(:,:,:)
        
        integer*4 :: ncID,varID
        integer*4, dimension(nf90_max_var_dims) :: dimIDs
        integer*4 :: Dim01,Dim02,Dim03,Dim04
        integer*4 :: TimeInd
        character :: FLAG !!! FLAG == 'A':allocate, FLAG == 'N',noallocate
        
        call NC_Check(nf90_open(trim(FPath)//trim(FName),nf90_NoWrite,ncID));
        
        call NC_Check(nf90_inq_varid(ncID,trim(VarName),varID));
        
        call NC_Check(nf90_Inquire_Variable(ncID,varID,dimids = dimIDs))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(1), len=Dim01))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(2), len=Dim02))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(3), len=Dim03))
        
        IF( FLAG == 'A')THEN
        allocate(PutArray(Dim01,Dim02,Dim03));
        END IF
        
        call NC_Check( nf90_get_var(ncID,varID,PutArray,(/1,1,1,int4(TimeInd)/)) );
        
        call NC_Check(nf90_close(ncID))
        
        write(*,*) 'Reading '//trim(VarName)//' from '//trim(FPath)//trim(FName)
        write(*,*),''
        
    END SUBROUTINE NC_Read_4Dto3D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Read_5Dto4D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Read_5Dto4D(FPath,FName,VarName,PutArray,TimeInd,FLAG)
    
        character(len=*) :: FPath,FName,VarName
        real*8,allocatable :: PutArray(:,:,:,:)
        
        integer*4 :: ncID,varID
        integer*4, dimension(nf90_max_var_dims) :: dimIDs
        integer*4 :: Dim01,Dim02,Dim03,Dim04,Dim05
        integer*4 :: TimeInd
        character :: FLAG !!! FLAG == 'A':allocate, FLAG == 'N',noallocate
        
        call NC_Check(nf90_open(trim(FPath)//trim(FName),nf90_NoWrite,ncID));
        
        call NC_Check(nf90_inq_varid(ncID,trim(VarName),varID));
        
        call NC_Check(nf90_Inquire_Variable(ncID,varID,dimids = dimIDs))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(1), len=Dim01))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(2), len=Dim02))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(3), len=Dim03))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(4), len=Dim04))
        
        IF( FLAG == 'A')THEN
        allocate(PutArray(Dim01,Dim02,Dim03,Dim04));
        END IF
        
        call NC_Check( nf90_get_var(ncID,varID,PutArray,(/1,1,1,1,int4(TimeInd)/)) );
        
        call NC_Check(nf90_close(ncID))
        
        write(*,*) 'Reading '//trim(VarName)//' from '//trim(FPath)//trim(FName)
        write(*,*),''
        
    END SUBROUTINE NC_Read_5Dto4D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!
    !!! NC_Create_File
    !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Create_1D_File (Space 0D)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Create_1D_File &
               ( FPath,FName, Dim01_Val,Dim01_Name, Time_Unit, LeapYear_FLAG)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
       character(len=*) :: FPath,FName
       character(len=*) ::   Dim01_Name, Time_Unit  !!! dim unit
       real*8,allocatable :: Dim01_Val(:) !!! dim value
       integer*8 ::  Dim01_Num !!! dim num
       character ::  GridType
       
       integer*8,optional,intent(in) :: LeapYear_FLAG 
       !!! use Leap Year (1) or not (0) 
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       integer*4 :: ncFileID, Dim01_ID
       integer*4 :: Dim01_VarID
       
       ! get Dim Num
       Dim01_Num = size(Dim01_Val)
       
       ! Create the file
       call NC_Check(nf90_create(path = trim(FPath)//trim(FName), cmode = nf90_64bit_offset, ncid = ncFileID))
       
       ! Define the dimensions
       call NC_Check(nf90_def_dim(ncid = ncFileID, name = trim(Dim01_Name), len = nf90_unlimited, dimid = Dim01_ID))
       
       ! Define Dim Var
       call NC_Check(nf90_def_var(ncFileID, trim(Dim01_Name), nf90_double, dimids = Dim01_ID, varID = Dim01_VarID))
       call NC_Check(nf90_put_att(ncFileID, Dim01_VarID, 'units', trim(Time_Unit)))
       
       ! default use 365_day calendar 
       IF(present(LeapYear_FLAG))THEN
          IF(LeapYear_FLAG == 1)THEN
            call NC_Check(nf90_put_att(ncFileID, Dim01_VarID, 'calendar','gregorian'))
          ELSE
            call NC_Check(nf90_put_att(ncFileID, Dim01_VarID, 'calendar','365_day'))
          ENDIF
       ELSE
            call NC_Check(nf90_put_att(ncFileID, Dim01_VarID, 'calendar','365_day'))
       ENDIF
       
       ! end def
       call NC_Check(nf90_enddef(ncfileID))
       
       ! put Dim Var Value
       call NC_Check(nf90_put_var(ncFileID, Dim01_VarID, Dim01_Val))
       
       ! close file
       call NC_Check(nf90_close(ncfileID))
       
       write(*,*) 'Creating '//trim(FPath)//trim(FName)
       write(*,*),''
    
    END SUBROUTINE NC_Create_1D_File
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Create_2D_File (Space 1D)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Create_2D_File &
               ( FPath,FName, Dim01_Val, Dim02_Val, Dim01_Name, Dim02_Name,Time_Unit,LeapYear_FLAG)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
       character(len=*) :: FPath,FName
       character(len=*) ::   Dim01_Name, Dim02_Name, Time_Unit  !!! dim unit
       real*8,allocatable :: Dim01_Val(:),Dim02_Val(:) !!! dim value
       integer*8 ::  Dim01_Num,Dim02_Num !!! dim num
       character ::  GridType
       
       integer*8,optional,intent(in) :: LeapYear_FLAG 
       !!! use Leap Year (1) or not (0) 
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       integer*4 :: ncFileID, Dim01_ID, Dim02_ID
       integer*4 :: Dim01_VarID, Dim02_VarID
       
       ! get Dim Num
       Dim01_Num = size(Dim01_Val)
       Dim02_Num = size(Dim02_Val)
       
       ! Create the file
       call NC_Check(nf90_create(path = trim(FPath)//trim(FName), cmode = nf90_64bit_offset, ncid = ncFileID))
       
       ! Define the dimensions
       call NC_Check(nf90_def_dim(ncid = ncFileID, name = trim(Dim01_Name), len = int4(Dim01_Num), dimid = Dim01_ID))
       call NC_Check(nf90_def_dim(ncid = ncFileID, name = trim(Dim02_Name), len = nf90_unlimited, dimid = Dim02_ID))
       
       ! Define Dim Var
       call NC_Check(nf90_def_var(ncFileID, trim(Dim01_Name), nf90_double, dimids = Dim01_ID, varID = Dim01_VarID))
       call NC_Check(nf90_def_var(ncFileID, trim(Dim02_Name), nf90_double, dimids = Dim02_ID, varID = Dim02_VarID))
       call NC_Check(nf90_put_att(ncFileID, Dim02_VarID, 'units', trim(Time_Unit)))
       
       
       ! default use 365_day calendar 
       IF(present(LeapYear_FLAG))THEN
          IF(LeapYear_FLAG == 1)THEN
            call NC_Check(nf90_put_att(ncFileID, Dim02_VarID, 'calendar','gregorian'))
          ELSE
            call NC_Check(nf90_put_att(ncFileID, Dim02_VarID, 'calendar','365_day'))
          ENDIF
       ELSE
            call NC_Check(nf90_put_att(ncFileID, Dim02_VarID, 'calendar','365_day'))
       ENDIF
       
       ! end def
       call NC_Check(nf90_enddef(ncfileID))
       
       ! put Dim Var Value
       call NC_Check(nf90_put_var(ncFileID, Dim01_VarID, Dim01_Val))
       call NC_Check(nf90_put_var(ncFileID, Dim02_VarID, Dim02_Val))
       
       ! close file
       call NC_Check(nf90_close(ncfileID))
    
       write(*,*) 'Creating '//trim(FPath)//trim(FName)
       write(*,*),''
    
    END SUBROUTINE NC_Create_2D_File
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Create_3D_File (Space 2D)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Create_3D_File &
               ( FPath,FName, &
                 Dim01_Val, Dim02_Val, Dim03_Val, &
                 Dim01_Name, Dim02_Name, Dim03_Name, Time_Unit, GridType, LeapYear_FLAG)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
       character(len=*) :: FPath,FName
       character(len=*) ::   Dim01_Name, Dim02_Name, Dim03_Name, Time_Unit  !!! dim unit
       real*8,allocatable :: Dim01_Val(:),Dim02_Val(:),Dim03_Val(:) !!! dim value
       integer*8 ::  Dim01_Num,Dim02_Num,Dim03_Num !!! dim num
       character ::  GridType
       
       integer*8,optional,intent(in) :: LeapYear_FLAG 
       !!! use Leap Year (1) or not (0) 
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       integer*4 :: ncFileID, Dim01_ID, Dim02_ID, Dim03_ID
       integer*4 :: Dim01_VarID, Dim02_VarID, Dim03_VarID
       
       ! get Dim Num
       Dim01_Num = size(Dim01_Val)
       Dim02_Num = size(Dim02_Val)
       Dim03_Num = size(Dim03_Val)
       
       ! Create the file
       call NC_Check(nf90_create(path = trim(FPath)//trim(FName), cmode = nf90_64bit_offset, ncid = ncFileID))
       
       ! Define the dimensions
       call NC_Check(nf90_def_dim(ncid = ncFileID, name = trim(Dim01_Name), len = int4(Dim01_Num), dimid = Dim01_ID))
       call NC_Check(nf90_def_dim(ncid = ncFileID, name = trim(Dim02_Name), len = int4(Dim02_Num), dimid = Dim02_ID))
       call NC_Check(nf90_def_dim(ncid = ncFileID, name = trim(Dim03_Name), len = nf90_unlimited, dimid = Dim03_ID))
       
       ! Define Dim Var
       call NC_Check(nf90_def_var(ncFileID, trim(Dim01_Name), nf90_double, dimids = Dim01_ID, varID = Dim01_VarID))
       call NC_Check(nf90_def_var(ncFileID, trim(Dim02_Name), nf90_double, dimids = Dim02_ID, varID = Dim02_VarID))
       call NC_Check(nf90_def_var(ncFileID, trim(Dim03_Name), nf90_double, dimids = Dim03_ID, varID = Dim03_VarID))
       call NC_Check(nf90_put_att(ncFileID, Dim03_VarID, 'units', trim(Time_Unit)))
       
       
       ! default use 365_day calendar 
       IF(present(LeapYear_FLAG))THEN
          IF(LeapYear_FLAG == 1)THEN
            call NC_Check(nf90_put_att(ncFileID, Dim03_VarID, 'calendar','gregorian'))
          ELSE
            call NC_Check(nf90_put_att(ncFileID, Dim03_VarID, 'calendar','365_day'))
          ENDIF
       ELSE
            call NC_Check(nf90_put_att(ncFileID, Dim03_VarID, 'calendar','365_day'))
       ENDIF
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF(GridType == 'Y')THEN
           
          call NC_Check(nf90_put_att(ncFileID, Dim01_VarID, 'units', 'degrees_north'))
          call NC_Check(nf90_put_att(ncFileID, Dim01_VarID, 'long_name', 'Latitude'))
          
          call NC_Check(nf90_put_att(ncFileID, Dim02_VarID, 'units', 'degrees_east'))
          call NC_Check(nf90_put_att(ncFileID, Dim02_VarID, 'long_name', 'Longitude'))
          
       ELSEIF(GridType == 'X')THEN
           
          call NC_Check(nf90_put_att(ncFileID, Dim02_VarID, 'units', 'degrees_north'))
          call NC_Check(nf90_put_att(ncFileID, Dim02_VarID, 'long_name', 'Latitude'))
          
          call NC_Check(nf90_put_att(ncFileID, Dim01_VarID, 'units', 'degrees_east'))
          call NC_Check(nf90_put_att(ncFileID, Dim01_VarID, 'long_name', 'Longitude'))
           
       END IF
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       ! end def
       call NC_Check(nf90_enddef(ncfileID))
       
       ! put Dim Var Value
       call NC_Check(nf90_put_var(ncFileID, Dim01_VarID, Dim01_Val))
       call NC_Check(nf90_put_var(ncFileID, Dim02_VarID, Dim02_Val))
       call NC_Check(nf90_put_var(ncFileID, Dim03_VarID, Dim03_Val))
       
       ! close file
       call NC_Check(nf90_close(ncfileID))
       
       write(*,*) 'Creating '//trim(FPath)//trim(FName)
       write(*,*),''
    
    END SUBROUTINE NC_Create_3D_File
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Create_4D_File (Space 3D)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Create_4D_File &
               ( FPath,FName, &
                 Dim01_Val, Dim02_Val, Dim03_Val, Dim04_Val, &
                 Dim01_Name, Dim02_Name, Dim03_Name, Dim04_Name, Time_Unit, GridType, LeapYear_FLAG)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
       character(len=*) :: FPath,FName
       character(len=*) ::   Dim01_Name, Dim02_Name, Dim03_Name, Dim04_Name, Time_Unit  !!! dim unit
       real*8,allocatable :: Dim01_Val(:),Dim02_Val(:),Dim03_Val(:),Dim04_Val(:) !!! dim value
       integer*8 ::  Dim01_Num,Dim02_Num,Dim03_Num,Dim04_Num  !!! dim num
       character ::  GridType
       
       integer*8,optional,intent(in) :: LeapYear_FLAG 
       !!! use Leap Year (1) or not (0) 
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       integer*4 :: ncFileID, Dim01_ID, Dim02_ID, Dim03_ID, Dim04_ID
       integer*4 :: Dim01_VarID, Dim02_VarID, Dim03_VarID, Dim04_VarID
       
       ! get Dim Num
       Dim01_Num = size(Dim01_Val)
       Dim02_Num = size(Dim02_Val)
       Dim03_Num = size(Dim03_Val)
       Dim04_Num = size(Dim04_Val)
       
       ! Create the file
       call NC_Check(nf90_create(path = trim(FPath)//trim(FName), cmode = nf90_64bit_offset, ncid = ncFileID))
       
       ! Define the dimensions
       call NC_Check(nf90_def_dim(ncid = ncFileID, name = trim(Dim01_Name), len = int4(Dim01_Num), dimid = Dim01_ID))
       call NC_Check(nf90_def_dim(ncid = ncFileID, name = trim(Dim02_Name), len = int4(Dim02_Num), dimid = Dim02_ID))
       call NC_Check(nf90_def_dim(ncid = ncFileID, name = trim(Dim03_Name), len = int4(Dim03_Num), dimid = Dim03_ID))
       call NC_Check(nf90_def_dim(ncid = ncFileID, name = trim(Dim04_Name), len = nf90_unlimited, dimid = Dim04_ID))
       
       ! Define Dim Var
       call NC_Check(nf90_def_var(ncFileID, trim(Dim01_Name), nf90_double, dimids = Dim01_ID, varID = Dim01_VarID))
       call NC_Check(nf90_def_var(ncFileID, trim(Dim02_Name), nf90_double, dimids = Dim02_ID, varID = Dim02_VarID))
       call NC_Check(nf90_def_var(ncFileID, trim(Dim03_Name), nf90_double, dimids = Dim03_ID, varID = Dim03_VarID))
       call NC_Check(nf90_def_var(ncFileID, trim(Dim04_Name), nf90_double, dimids = Dim04_ID, varID = Dim04_VarID))
       call NC_Check(nf90_put_att(ncFileID, Dim04_VarID, 'units', trim(Time_Unit)))
       
       
       ! default use 365_day calendar
       IF(present(LeapYear_FLAG))THEN
          IF(LeapYear_FLAG == 1)THEN
            call NC_Check(nf90_put_att(ncFileID, Dim04_VarID, 'calendar','gregorian'))
          ELSE
            call NC_Check(nf90_put_att(ncFileID, Dim04_VarID, 'calendar','365_day'))
          ENDIF
       ELSE
            call NC_Check(nf90_put_att(ncFileID, Dim04_VarID, 'calendar','365_day'))
       ENDIF
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF(GridType == 'Y')THEN
           
          call NC_Check(nf90_put_att(ncFileID, Dim01_VarID, 'units', 'degrees_north'))
          call NC_Check(nf90_put_att(ncFileID, Dim01_VarID, 'long_name', 'Latitude'))
          
          call NC_Check(nf90_put_att(ncFileID, Dim02_VarID, 'units', 'degrees_east'))
          call NC_Check(nf90_put_att(ncFileID, Dim02_VarID, 'long_name', 'Longitude'))
          
       ELSEIF(GridType == 'X')THEN
           
          call NC_Check(nf90_put_att(ncFileID, Dim02_VarID, 'units', 'degrees_north'))
          call NC_Check(nf90_put_att(ncFileID, Dim02_VarID, 'long_name', 'Latitude'))
          
          call NC_Check(nf90_put_att(ncFileID, Dim01_VarID, 'units', 'degrees_east'))
          call NC_Check(nf90_put_att(ncFileID, Dim01_VarID, 'long_name', 'Longitude'))
           
       END IF
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       ! end def
       call NC_Check(nf90_enddef(ncfileID))
       
       ! put Dim Var Value
       call NC_Check(nf90_put_var(ncFileID, Dim01_VarID, Dim01_Val))
       call NC_Check(nf90_put_var(ncFileID, Dim02_VarID, Dim02_Val))
       call NC_Check(nf90_put_var(ncFileID, Dim03_VarID, Dim03_Val))
       call NC_Check(nf90_put_var(ncFileID, Dim04_VarID, Dim04_Val))
       
       ! close file
       call NC_Check(nf90_close(ncfileID))
       
       write(*,*) 'Creating '//trim(FPath)//trim(FName)
       write(*,*),''
    
    END SUBROUTINE NC_Create_4D_File
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Create_5D_File (Space 4D)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Create_5D_File &
              (FPath,FName, &
               Dim01_Val, Dim02_Val, Dim03_Val, Dim04_Val, Dim05_Val, &
               Dim01_Name, Dim02_Name, Dim03_Name, Dim04_Name, Dim05_Name, Time_Unit, GridType, LeapYear_FLAG)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
       character(len=*) :: FPath,FName
       character(len=*) ::   Dim01_Name, Dim02_Name, Dim03_Name, Dim04_Name, Dim05_Name, Time_Unit  !!! dim unit
       real*8,allocatable :: Dim01_Val(:),Dim02_Val(:),Dim03_Val(:),Dim04_Val(:), Dim05_Val(:) !!! dim value
       integer*8 ::  Dim01_Num,Dim02_Num,Dim03_Num,Dim04_Num,Dim05_Num  !!! dim num
       character ::  GridType
       
       integer*8,optional,intent(in) :: LeapYear_FLAG 
       !!! use Leap Year (1) or not (0) 
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       integer*4 :: ncFileID, Dim01_ID, Dim02_ID, Dim03_ID, Dim04_ID, Dim05_ID
       integer*4 :: Dim01_VarID, Dim02_VarID, Dim03_VarID, Dim04_VarID, Dim05_VarID
       
       ! get Dim Num
       Dim01_Num = size(Dim01_Val)
       Dim02_Num = size(Dim02_Val)
       Dim03_Num = size(Dim03_Val)
       Dim04_Num = size(Dim04_Val)
       Dim05_Num = size(Dim05_Val)
       
       ! Create the file
       call NC_Check(nf90_create(path = trim(FPath)//trim(FName), cmode = nf90_64bit_offset, ncid = ncFileID))
       
       ! Define the dimensions
       call NC_Check(nf90_def_dim(ncid = ncFileID, name = trim(Dim01_Name), len = int4(Dim01_Num), dimid = Dim01_ID))
       call NC_Check(nf90_def_dim(ncid = ncFileID, name = trim(Dim02_Name), len = int4(Dim02_Num), dimid = Dim02_ID))
       call NC_Check(nf90_def_dim(ncid = ncFileID, name = trim(Dim03_Name), len = int4(Dim03_Num), dimid = Dim03_ID))
       call NC_Check(nf90_def_dim(ncid = ncFileID, name = trim(Dim04_Name), len = int4(Dim04_Num), dimid = Dim04_ID))
       call NC_Check(nf90_def_dim(ncid = ncFileID, name = trim(Dim05_Name), len = nf90_unlimited,  dimid = Dim05_ID))
       
       ! Define Dim Var
       call NC_Check(nf90_def_var(ncFileID, trim(Dim01_Name), nf90_double, dimids = Dim01_ID, varID = Dim01_VarID))
       call NC_Check(nf90_def_var(ncFileID, trim(Dim02_Name), nf90_double, dimids = Dim02_ID, varID = Dim02_VarID))
       call NC_Check(nf90_def_var(ncFileID, trim(Dim03_Name), nf90_double, dimids = Dim03_ID, varID = Dim03_VarID))
       call NC_Check(nf90_def_var(ncFileID, trim(Dim04_Name), nf90_double, dimids = Dim04_ID, varID = Dim04_VarID))
       call NC_Check(nf90_def_var(ncFileID, trim(Dim05_Name), nf90_double, dimids = Dim05_ID, varID = Dim05_VarID))
       call NC_Check(nf90_put_att(ncFileID, Dim05_VarID, 'units', trim(Time_Unit)))
       
       
       ! default use 365_day calendar
       IF(present(LeapYear_FLAG))THEN
          IF(LeapYear_FLAG == 1)THEN
            call NC_Check(nf90_put_att(ncFileID, Dim05_VarID, 'calendar','gregorian'))
          ELSE
            call NC_Check(nf90_put_att(ncFileID, Dim05_VarID, 'calendar','365_day'))
          ENDIF
       ELSE
            call NC_Check(nf90_put_att(ncFileID, Dim05_VarID, 'calendar','365_day'))
       ENDIF
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF(GridType == 'Y')THEN
           
          call NC_Check(nf90_put_att(ncFileID, Dim01_VarID, 'units', 'degrees_north'))
          call NC_Check(nf90_put_att(ncFileID, Dim01_VarID, 'long_name', 'Latitude'))
          
          call NC_Check(nf90_put_att(ncFileID, Dim02_VarID, 'units', 'degrees_east'))
          call NC_Check(nf90_put_att(ncFileID, Dim02_VarID, 'long_name', 'Longitude'))
          
       ELSEIF(GridType == 'X')THEN
           
          call NC_Check(nf90_put_att(ncFileID, Dim02_VarID, 'units', 'degrees_north'))
          call NC_Check(nf90_put_att(ncFileID, Dim02_VarID, 'long_name', 'Latitude'))
          
          call NC_Check(nf90_put_att(ncFileID, Dim01_VarID, 'units', 'degrees_east'))
          call NC_Check(nf90_put_att(ncFileID, Dim01_VarID, 'long_name', 'Longitude'))
           
       END IF
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       ! end def
       call NC_Check(nf90_enddef(ncfileID))
       
       ! put Dim Var Value
       call NC_Check(nf90_put_var(ncFileID, Dim01_VarID, Dim01_Val))
       call NC_Check(nf90_put_var(ncFileID, Dim02_VarID, Dim02_Val))
       call NC_Check(nf90_put_var(ncFileID, Dim03_VarID, Dim03_Val))
       call NC_Check(nf90_put_var(ncFileID, Dim04_VarID, Dim04_Val))
       call NC_Check(nf90_put_var(ncFileID, Dim05_VarID, Dim05_Val))
       
       ! close file
       call NC_Check(nf90_close(ncfileID))
       
       write(*,*) 'Creating '//trim(FPath)//trim(FName)
       write(*,*),''
    
    END SUBROUTINE NC_Create_5D_File
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!
    !!! NC_Def_Var
    !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Def_1D_Var
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Def_1D_Var(FPath,FName,VarName,VarUnit, Dim01_Name, Double_FLAG)
    
    
       character(len=*) :: FPath,FName
       character(len=*) :: VarName, VarUnit
       character(len=*) :: Dim01_Name
       
       integer*4 :: ncFileID, Dim01_ID, VarID
       integer*4, optional, intent(in) :: Double_FLAG !!! double output (1) or not (0)
       
       !open and redef
       call NC_Check(nf90_open(path = trim(FPath)//trim(FName), mode = nf90_write, ncid = ncFileID))
       call NC_Check(nf90_redef(ncFileID))
       
       !inquire dimid
       call NC_Check(nf90_inq_dimid(ncFileID,trim(Dim01_Name),Dim01_ID))
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! default single precision (float) output
       IF( present(Double_FLAG) ) THEN
           IF( Double_FLAG == 1) THEN
           call NC_Check(nf90_def_var(ncid = ncFileID, name = trim(VarName), xtype = nf90_double, &
                                      dimids = Dim01_ID, varID = VarID) )
           ELSE
           call NC_Check(nf90_def_var(ncid = ncFileID, name = trim(VarName), xtype = nf90_float, &
                                      dimids = Dim01_ID, varID = VarID) )
           END IF
       ELSE
           call NC_Check(nf90_def_var(ncid = ncFileID, name = trim(VarName), xtype = nf90_float, &
                                      dimids = Dim01_ID, varID = VarID) )
       END IF
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       call NC_Check(nf90_put_att(ncFileID, VarID, "units",trim(VarUnit)) )
       
       !end def and close
       call NC_Check(nf90_enddef(ncFileID))
       call NC_Check(nf90_close(ncFileID))
    
    
    END SUBROUTINE NC_Def_1D_Var
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                         
    !!! NC_Def_2D_Var
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Def_2D_Var(FPath,FName,VarName,VarUnit, Dim01_Name, Dim02_Name, Double_FLAG)
    
    
       character(len=*) :: FPath,FName
       character(len=*) :: VarName, VarUnit
       character(len=*) :: Dim01_Name,Dim02_Name
       
       integer*4 :: ncFileID, Dim01_ID, Dim02_ID, VarID
       integer*4, optional, intent(in) :: Double_FLAG !!! double output (1) or not (0)
       
       !open and redef
       call NC_Check(nf90_open(path = trim(FPath)//trim(FName), mode = nf90_write, ncid = ncFileID))
       call NC_Check(nf90_redef(ncFileID))
       
       !inquire dimid
       call NC_Check(nf90_inq_dimid(ncFileID,trim(Dim01_Name),Dim01_ID))
       call NC_Check(nf90_inq_dimid(ncFileID,trim(Dim02_Name),Dim02_ID))
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! default single precision (float) output
       IF( present(Double_FLAG) ) THEN
           IF( Double_FLAG == 1) THEN
           call NC_Check(nf90_def_var(ncid = ncFileID, name = trim(VarName), xtype = nf90_double, &
                                      dimids = (/Dim01_ID, Dim02_ID /), varID = VarID) )
           ELSE
           call NC_Check(nf90_def_var(ncid = ncFileID, name = trim(VarName), xtype = nf90_float, &
                                      dimids = (/Dim01_ID, Dim02_ID /), varID = VarID) )
           END IF
       ELSE
           call NC_Check(nf90_def_var(ncid = ncFileID, name = trim(VarName), xtype = nf90_float, &
                                      dimids = (/Dim01_ID, Dim02_ID /), varID = VarID) )
       END IF
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       call NC_Check(nf90_put_att(ncFileID, VarID, "units",trim(VarUnit)) )
       
       !end def and close
       call NC_Check(nf90_enddef(ncFileID))
       call NC_Check(nf90_close(ncFileID))
    
    
    END SUBROUTINE NC_Def_2D_Var
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Def_3D_Var
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Def_3D_Var(FPath,FName,VarName,VarUnit, &
                             Dim01_Name, Dim02_Name, Dim03_Name, Double_FLAG)
    
       character(len=*) :: FPath,FName
       character(len=*) :: VarName, VarUnit
       character(len=*) :: Dim01_Name,Dim02_Name,Dim03_Name !!! dim name: usually (lat,lon,time)
       
       integer*4 :: ncFileID, Dim01_ID, Dim02_ID, Dim03_ID, VarID
       integer*4, optional, intent(in) :: Double_FLAG !!! double output (1) or not (0)
       
       !open and redef
       call NC_Check(nf90_open(path = trim(FPath)//trim(FName), mode = nf90_write, ncid = ncFileID))
       call NC_Check(nf90_redef(ncFileID))
       
       !inquire dimid
       call NC_Check(nf90_inq_dimid(ncFileID,trim(Dim01_Name),Dim01_ID))
       call NC_Check(nf90_inq_dimid(ncFileID,trim(Dim02_Name),Dim02_ID))
       call NC_Check(nf90_inq_dimid(ncFileID,trim(Dim03_Name),Dim03_ID))
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! default single precision (float) output
       IF( present(Double_FLAG) ) THEN
           IF( Double_FLAG == 1) THEN
           call NC_Check(nf90_def_var(ncid = ncFileID, name = trim(VarName), xtype = nf90_double, &
                                      dimids = (/Dim01_ID, Dim02_ID, Dim03_ID /), varID = VarID) )
           ELSE
           call NC_Check(nf90_def_var(ncid = ncFileID, name = trim(VarName), xtype = nf90_float, &
                                      dimids = (/Dim01_ID, Dim02_ID, Dim03_ID /), varID = VarID) )
           END IF
       ELSE
           call NC_Check(nf90_def_var(ncid = ncFileID, name = trim(VarName), xtype = nf90_float, &
                                      dimids = (/Dim01_ID, Dim02_ID, Dim03_ID /), varID = VarID) )
       END IF
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       call NC_Check(nf90_put_att(ncFileID, VarID, "units",trim(VarUnit)) )
       
       !end def and close
       call NC_Check(nf90_enddef(ncFileID))
       call NC_Check(nf90_close(ncFileID))
    
    
    END SUBROUTINE NC_Def_3D_Var
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
    !!! NC_Def_4D_Var
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Def_4D_Var(FPath,FName,VarName,VarUnit, &
                             Dim01_Name,Dim02_Name,Dim03_Name, Dim04_Name, Double_FLAG)
    
    
       character(len=*) :: FPath,FName
       character(len=*) :: VarName, VarUnit
       character(len=*) :: Dim01_Name,Dim02_Name,Dim03_Name,Dim04_Name 
       !!! dim name: usually (lat,lon,depth,time)
       
       integer*4 :: ncFileID, Dim01_ID, Dim02_ID, Dim03_ID, Dim04_ID, VarID
       integer*4, optional, intent(in) :: Double_FLAG !!! double output (1) or not (0)
       
       !open and redef
       call NC_Check(nf90_open(path = trim(FPath)//trim(FName), mode = nf90_write, ncid = ncFileID))
       call NC_Check(nf90_redef(ncFileID))
       
       !inquire dimid
       call NC_Check(nf90_inq_dimid(ncFileID,trim(Dim01_Name),Dim01_ID))
       call NC_Check(nf90_inq_dimid(ncFileID,trim(Dim02_Name),Dim02_ID))
       call NC_Check(nf90_inq_dimid(ncFileID,trim(Dim03_Name),Dim03_ID))
       call NC_Check(nf90_inq_dimid(ncFileID,trim(Dim04_Name),Dim04_ID))
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! default single precision (float) output
       IF( present(Double_FLAG) ) THEN
           IF( Double_FLAG == 1) THEN
           call NC_Check(nf90_def_var(ncid = ncFileID, name = trim(VarName), xtype = nf90_double, &
                         dimids = (/Dim01_ID, Dim02_ID, Dim03_ID, Dim04_ID /), varID = VarID) )
           ELSE
           call NC_Check(nf90_def_var(ncid = ncFileID, name = trim(VarName), xtype = nf90_float, &
                         dimids = (/Dim01_ID, Dim02_ID, Dim03_ID, Dim04_ID /), varID = VarID) )
           END IF
       ELSE
           call NC_Check(nf90_def_var(ncid = ncFileID, name = trim(VarName), xtype = nf90_float, &
                         dimids = (/Dim01_ID, Dim02_ID, Dim03_ID, Dim04_ID /), varID = VarID) )
       END IF
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       
       call NC_Check(nf90_put_att(ncFileID, VarID, "units",trim(VarUnit)) )
       
       !end def and close
       call NC_Check(nf90_enddef(ncFileID))
       call NC_Check(nf90_close(ncFileID))
    
    
    END SUBROUTINE NC_Def_4D_Var
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
    !!! NC_Def_5D_Var
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Def_5D_Var(FPath,FName,VarName,VarUnit, &
                             Dim01_Name,Dim02_Name,Dim03_Name, Dim04_Name, Dim05_Name, Double_FLAG)
    
    
       character(len=*) :: FPath,FName
       character(len=*) :: VarName, VarUnit
       character(len=*) :: Dim01_Name,Dim02_Name,Dim03_Name,Dim04_Name, Dim05_Name
       !!! dim name: usually (lat,lon,depth,time)
       
       integer*4 :: ncFileID, Dim01_ID, Dim02_ID, Dim03_ID, Dim04_ID, Dim05_ID, VarID
       integer*4, optional, intent(in) :: Double_FLAG !!! double output (1) or not (0)
       
       !open and redef
       call NC_Check(nf90_open(path = trim(FPath)//trim(FName), mode = nf90_write, ncid = ncFileID))
       call NC_Check(nf90_redef(ncFileID))
       
       !inquire dimid
       call NC_Check(nf90_inq_dimid(ncFileID,trim(Dim01_Name),Dim01_ID))
       call NC_Check(nf90_inq_dimid(ncFileID,trim(Dim02_Name),Dim02_ID))
       call NC_Check(nf90_inq_dimid(ncFileID,trim(Dim03_Name),Dim03_ID))
       call NC_Check(nf90_inq_dimid(ncFileID,trim(Dim04_Name),Dim04_ID))
       call NC_Check(nf90_inq_dimid(ncFileID,trim(Dim05_Name),Dim04_ID))
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! default single precision (float) output
       IF( present(Double_FLAG) ) THEN
           IF( Double_FLAG == 1) THEN
           call NC_Check(nf90_def_var(ncid = ncFileID, name = trim(VarName), xtype = nf90_double, &
                         dimids = (/Dim01_ID, Dim02_ID, Dim03_ID, Dim04_ID, Dim05_ID /), varID = VarID) )
           ELSE
           call NC_Check(nf90_def_var(ncid = ncFileID, name = trim(VarName), xtype = nf90_float, &
                         dimids = (/Dim01_ID, Dim02_ID, Dim03_ID, Dim04_ID, Dim05_ID /), varID = VarID) )
           END IF
       ELSE
           call NC_Check(nf90_def_var(ncid = ncFileID, name = trim(VarName), xtype = nf90_float, &
                         dimids = (/Dim01_ID, Dim02_ID, Dim03_ID, Dim04_ID, Dim05_ID /), varID = VarID) )
       END IF
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       
       call NC_Check(nf90_put_att(ncFileID, VarID, "units",trim(VarUnit)) )
       
       !end def and close
       call NC_Check(nf90_enddef(ncFileID))
       call NC_Check(nf90_close(ncFileID))
    
    
    END SUBROUTINE NC_Def_5D_Var
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!
    !!! NC_Write_Mul_IT
    !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Write_1Dto1D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Write_1Dto1D(FPath, FName,  VarName, Value, StartIT, CountT)
    
        character(len=*) :: FPath,FName,VarName
        integer*4 ::  StartIT, CountT
        real*8,allocatable :: Value(:)
        
        integer*4 :: ncID, rhVarId
        integer*4, dimension(nf90_max_var_dims) :: dimIDs
        integer*4 :: Dim01
        
        call NC_Check(nf90_open(path = trim(FPath)//trim(FName), mode = nf90_write, ncid = ncID))
        call NC_Check(nf90_inq_varid(ncid, trim(VarName), rhVarId))
        
        ! get space Dim
        call NC_Check(nf90_Inquire_Variable(ncID, rhVarId, dimids = dimIDs))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(1), len=Dim01))
        
        ! put var
        IF(StartIT<=0)THEN
          call NC_Check(nf90_put_var(ncID,rhVarID,Value));
        ELSE
          call NC_Check(nf90_put_var(ncID,rhVarID,Value,&
                        (/int4(StartIT)/),(/int4(CountT)/)))
        END IF
        
        ! close
        call NC_Check(nf90_close(ncID))
        
        !write(*,*) 'Writing '//trim(VarName)//' to '//trim(FPath)//trim(FName)
        !write(*,*),''
        
    END SUBROUTINE NC_Write_1Dto1D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Write_2Dto2D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Write_2Dto2D(FPath, FName,VarName, Value, StartIT, CountT)
    
        character(len=*) :: FPath,FName,VarName
        integer*4 ::  StartIT, CountT
        real*8,allocatable :: Value(:,:)
        
        integer*4 :: ncID, rhVarId
        integer*4, dimension(nf90_max_var_dims) :: dimIDs
        integer*4 :: Dim01,Dim02
        
        call NC_Check(nf90_open(path = trim(FPath)//trim(FName), mode = nf90_write, ncid = ncID))
        call NC_Check(nf90_inq_varid(ncid, trim(VarName), rhVarId))
        
        ! get space Dim
        call NC_Check(nf90_Inquire_Variable(ncID, rhVarId, dimids = dimIDs))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(1), len=Dim01))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(2), len=Dim02))
        
        ! put var
        IF(StartIT<=0)THEN
          call NC_Check(nf90_put_var(ncID,rhVarID,Value));
        ELSE
          call NC_Check(nf90_put_var(ncID,rhVarID,Value,&
                        (/1,int4(StartIT)/),(/Dim01,int4(CountT)/)))
        END IF
        
        ! close
        call NC_Check(nf90_close(ncID))
        
        !write(*,*) 'Writing '//trim(VarName)//' to '//trim(FPath)//trim(FName)
        !write(*,*),''
        
    END SUBROUTINE NC_Write_2Dto2D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Write_3Dto3D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Write_3Dto3D(FPath, FName,VarName, Value, StartIT, CountT)
    
        character(len=*) :: FPath,FName,VarName
        integer*4 ::  StartIT, CountT
        real*8,allocatable :: Value(:,:,:)
        
        integer*4 :: ncID, rhVarId
        integer*4, dimension(nf90_max_var_dims) :: dimIDs
        integer*4 :: Dim01,Dim02,Dim03
        
        call NC_Check(nf90_open(path = trim(FPath)//trim(FName), mode = nf90_write, ncid = ncID))
        call NC_Check(nf90_inq_varid(ncid, trim(VarName), rhVarId))
        
        ! get space Dim
        call NC_Check(nf90_Inquire_Variable(ncID, rhVarId, dimids = dimIDs))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(1), len=Dim01))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(2), len=Dim02))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(3), len=Dim03))
        
        ! put var
        IF(StartIT<=0)THEN
          call NC_Check(nf90_put_var(ncID,rhVarID,Value));
        ELSE
          call NC_Check(nf90_put_var(ncID,rhVarID,Value,&
                        (/1,1,int4(StartIT)/),(/Dim01,Dim02,int4(CountT)/)))
        END IF
        
        ! close
        call NC_Check(nf90_close(ncID))
        
        !write(*,*) 'Writing '//trim(VarName)//' to '//trim(FPath)//trim(FName)
        !write(*,*),''
        
    END SUBROUTINE NC_Write_3Dto3D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Write_4Dto4D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Write_4Dto4D(FPath, FName,VarName, Value, StartIT, CountT)
    
        character(len=*) :: FPath,FName,VarName
        integer*4 :: StartIT,CountT
        real*8,allocatable :: Value(:,:,:,:)
        
        integer*4 :: ncID, rhVarId
        integer*4, dimension(nf90_max_var_dims) :: dimIDs
        integer*4 :: Dim01,Dim02,Dim03,Dim04
        
        call NC_Check(nf90_open(path = trim(FPath)//trim(FName), mode = nf90_write, ncid = ncID))
        call NC_Check(nf90_inq_varid(ncid, trim(VarName), rhVarId))
        
        ! get space Dim
        call NC_Check(nf90_Inquire_Variable(ncID, rhVarId, dimids = dimIDs))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(1), len=Dim01))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(2), len=Dim02))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(3), len=Dim03))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(4), len=Dim04))
        
        ! put var
        IF(StartIT<=0)THEN
          call NC_Check(nf90_put_var(ncID,rhVarID,Value));
        ELSE
          call NC_Check(nf90_put_var(ncID,rhVarID,Value,&
                        (/1,1,1,int4(StartIT)/),(/Dim01,Dim02,Dim03,int4(CountT)/)))
        END IF
        
        ! close
        call NC_Check(nf90_close(ncID))
        
        !write(*,*) 'Writing '//trim(VarName)//' to '//trim(FPath)//trim(FName)
        !write(*,*),''
        
    END SUBROUTINE NC_Write_4Dto4D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Write_5Dto5D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Write_5Dto5D(FPath, FName,VarName, Value, StartIT, CountT)
    
        character(len=*) :: FPath,FName,VarName
        integer*4 :: StartIT,CountT
        real*8,allocatable :: Value(:,:,:,:,:)
        
        integer*4 :: ncID, rhVarId
        integer*4, dimension(nf90_max_var_dims) :: dimIDs
        integer*4 :: Dim01,Dim02,Dim03,Dim04,Dim05
        
        call NC_Check(nf90_open(path = trim(FPath)//trim(FName), mode = nf90_write, ncid = ncID))
        call NC_Check(nf90_inq_varid(ncid, trim(VarName), rhVarId))
        
        ! get space Dim
        call NC_Check(nf90_Inquire_Variable(ncID, rhVarId, dimids = dimIDs))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(1), len=Dim01))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(2), len=Dim02))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(3), len=Dim03))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(4), len=Dim04))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(5), len=Dim05))
        
        ! put var
        IF(StartIT<=0)THEN
          call NC_Check(nf90_put_var(ncID,rhVarID,Value));
        ELSE
          call NC_Check(nf90_put_var(ncID,rhVarID,Value,&
                        (/1,1,1,1,int4(StartIT)/),(/Dim01,Dim02,Dim03,Dim04,int4(CountT)/)))
        END IF
        
        ! close
        call NC_Check(nf90_close(ncID))
        
        !write(*,*) 'Writing '//trim(VarName)//' to '//trim(FPath)//trim(FName)
        !write(*,*),''
        
    END SUBROUTINE NC_Write_5Dto5D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!
    !!! NC_Write_One_IT
    !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Write_0Dto1D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Write_0Dto1D(FPath, FName,  VarName,Value, TimeInd)
    
        character(len=*) :: FPath,FName,VarName
        integer*4 :: TimeInd
        real*8 :: Value
        
        integer*4 :: ncID, rhVarId
        integer*4, dimension(nf90_max_var_dims) :: dimIDs
        integer*4 :: Dim01
        
        call NC_Check(nf90_open(path = trim(FPath)//trim(FName), mode = nf90_write, ncid = ncID))
        call NC_Check(nf90_inq_varid(ncid, trim(VarName), rhVarId))
        
        ! get space Dim
        call NC_Check(nf90_Inquire_Variable(ncID, rhVarId, dimids = dimIDs))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(1), len=Dim01))
       
        ! put var
        call NC_Check(nf90_put_var(ncID,rhVarID,Value,(/int4(TimeInd)/)) );
        call NC_Check(nf90_close(ncID))
        
        !write(*,*) 'Writing '//trim(VarName)//' to '//trim(FPath)//trim(FName)
        !write(*,*),''
        
    END SUBROUTINE NC_Write_0Dto1D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Write_1Dto2D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Write_1Dto2D(FPath, FName, VarName, Value, TimeInd)
    
        character(len=*) :: FPath,FName,VarName
        integer*4 :: TimeInd
        real*8,allocatable :: Value(:)
        
        integer*4 :: ncID, rhVarId
        integer*4, dimension(nf90_max_var_dims) :: dimIDs
        integer*4 :: Dim01,Dim02
        
        call NC_Check(nf90_open(path = trim(FPath)//trim(FName), mode = nf90_write, ncid = ncID))
        call NC_Check(nf90_inq_varid(ncid, trim(VarName), rhVarId))
        
        ! get space Dim
        call NC_Check(nf90_Inquire_Variable(ncID, rhVarId, dimids = dimIDs))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(1), len=Dim01))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(2), len=Dim02))
        
        ! put var
        call NC_Check(nf90_put_var(ncID,rhVarID,Value,(/1,int4(TimeInd)/)) );
        call NC_Check(nf90_close(ncID))
        
        !write(*,*) 'Writing '//trim(VarName)//' to '//trim(FPath)//trim(FName)
        !write(*,*),''
        
    END SUBROUTINE NC_Write_1Dto2D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Write_2Dto3D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Write_2Dto3D(FPath, FName, VarName, Value, TimeInd)
    
        character(len=*) :: FPath,FName,VarName
        integer*4 :: TimeInd
        real*8,allocatable :: Value(:,:)
        
        integer*4 :: ncID, rhVarId
        integer*4, dimension(nf90_max_var_dims) :: dimIDs
        integer*4 :: Dim01,Dim02,Dim03
        
        call NC_Check(nf90_open(path = trim(FPath)//trim(FName), mode = nf90_write, ncid = ncID))
        call NC_Check(nf90_inq_varid(ncid, trim(VarName), rhVarId))
        
        ! get space Dim
        call NC_Check(nf90_Inquire_Variable(ncID, rhVarId, dimids = dimIDs))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(1), len=Dim01))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(2), len=Dim02))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(3), len=Dim03))
        
        ! put var
        call NC_Check(nf90_put_var(ncID,rhVarID,Value,(/1,1,int4(TimeInd)/)) );
        call NC_Check(nf90_close(ncID))
        
        !write(*,*) 'Writing '//trim(VarName)//' to '//trim(FPath)//trim(FName)
        !write(*,*),''
        
    END SUBROUTINE NC_Write_2Dto3D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Write_3Dto4D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Write_3Dto4D(FPath, FName,  VarName, Value,TimeInd)
    
        character(len=*) :: FPath,FName,VarName
        integer*4 :: TimeInd
        real*8,allocatable :: Value(:,:,:)
        
        integer*4 :: ncID, rhVarId
        integer*4, dimension(nf90_max_var_dims) :: dimIDs
        integer*4 :: Dim01,Dim02,Dim03,Dim04
        
        call NC_Check(nf90_open(path = trim(FPath)//trim(FName), mode = nf90_write, ncid = ncID))
        call NC_Check(nf90_inq_varid(ncid, trim(VarName), rhVarId))
        
        ! get space Dim
        call NC_Check(nf90_Inquire_Variable(ncID, rhVarId, dimids = dimIDs))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(1), len=Dim01))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(2), len=Dim02))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(3), len=Dim03))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(4), len=Dim04))
        
        ! put var
        call NC_Check(nf90_put_var(ncID,rhVarID,Value,(/1,1,1,int4(TimeInd)/)) );
        call NC_Check(nf90_close(ncID))
        
        !write(*,*) 'Writing '//trim(VarName)//' to '//trim(FPath)//trim(FName)
        !write(*,*),''
        
    END SUBROUTINE NC_Write_3Dto4D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! NC_Write_4Dto5D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE NC_Write_4Dto5D(FPath, FName,  VarName, Value,TimeInd)
    
        character(len=*) :: FPath,FName,VarName
        integer*4 :: TimeInd
        real*8,allocatable :: Value(:,:,:,:)
        
        integer*4 :: ncID, rhVarId
        integer*4, dimension(nf90_max_var_dims) :: dimIDs
        integer*4 :: Dim01,Dim02,Dim03,Dim04,Dim05
        
        call NC_Check(nf90_open(path = trim(FPath)//trim(FName), mode = nf90_write, ncid = ncID))
        call NC_Check(nf90_inq_varid(ncid, trim(VarName), rhVarId))
        
        ! get space Dim
        call NC_Check(nf90_Inquire_Variable(ncID, rhVarId, dimids = dimIDs))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(1), len=Dim01))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(2), len=Dim02))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(3), len=Dim03))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(4), len=Dim04))
        call NC_Check(nf90_Inquire_Dimension(ncID,dimIDs(5), len=Dim05))
        
        ! put var
        call NC_Check(nf90_put_var(ncID,rhVarID,Value,(/1,1,1,1,int4(TimeInd)/)) );
        call NC_Check(nf90_close(ncID))
        
        !write(*,*) 'Writing '//trim(VarName)//' to '//trim(FPath)//trim(FName)
        !write(*,*),''
        
    END SUBROUTINE NC_Write_4Dto5D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
END MODULE MO_NCTOOLS