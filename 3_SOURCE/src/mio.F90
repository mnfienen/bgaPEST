

module model_input_output

! -- Modules which are USED by this module are defined.

      use ERROR_MESSAGE
      use UTILITIES
      implicit none

      private

! -- The structure that is associated with this module is defined.

      type MIO_STRUC

        private

! -- Size

        integer                            :: numinfile=0  ! Number of input files
        integer                            :: numoutfile=0 ! Number of output files
        integer                            :: npar=0       ! Number of parameters
        integer                            :: nobs=0       ! Number of observations

! -- Counters and flags

        integer                            :: mcall                ! number of model calls
        integer                            :: template_status=0    ! status of template loading and reading MD???
        integer                            :: instruction_status=0 ! status of instruction loading and reading MD???

! -- Instruction-related variables

        integer                                                  :: asize=0      ! size of "a" 
        integer                                                  :: numl=0       ! size of "ll"
        integer                                                  :: ninstr=0     ! size of "lcins" 
        integer, dimension(:), pointer              :: ll         ! holds line advance data
        integer, dimension(:), pointer              :: lcins      ! pointer to instructions
        integer, dimension(:), pointer              :: obsn1(:)   ! stores part of observation instruction
        integer, dimension(:), pointer              :: obsn2(:)   ! stores part of observation instruction
        integer, dimension(:), pointer              :: iiobs(:)   ! stores part of observation instruction
        character (len=1), dimension(:), pointer    :: a          ! holds compressed instruction set
        character (len=20), dimension(:), pointer   :: aobs       ! name of each observation

! -- Template-related variables

        integer                                                  :: precis     ! precision protocol
        integer                                                  :: nopnt      ! decimal point protocol
        integer, dimension(:), pointer              :: nw         ! minimum word length of a parameter  MD???
        character (len=1),   dimension(:), pointer  :: mrkdel     ! marker delimiters                   MD???
        character (len=1),   dimension(:), pointer  :: pardel     ! parameter delimiters                MD???
        character (len=23),  dimension(:), pointer  :: pword      ! word for each parameter             MD???
        character (len=12),  dimension(:), pointer  :: apar       ! name of each parameter

! -- Filenames

        character (len=256)                         :: workdir=' '    ! model working directory
        character (len=256), dimension(:), pointer  :: tempfile       ! template files
        character (len=256), dimension(:), pointer  :: modinfile      ! model input files
        character (len=256), dimension(:), pointer  :: insfile        ! instruction files
        character (len=256), dimension(:), pointer  :: modoutfile     ! model output files

! -- Command to run model

        character (len=256)                :: comline        ! Command to run model

! -- Function call counters.

        integer                            :: ic_initialise=0
        integer                            :: ic_initialise_model_input=0
        integer                            :: ic_initialise_model_output=0
        integer                            :: ic_initialise_parameters=0
        integer                            :: ic_initialise_observations=0
        integer                            :: ic_set_number_precision=0
        integer                            :: ic_set_number_decpoint=0
        integer                            :: ic_set_model_command=0
        integer                            :: ic_get_model_command=0
        integer                            :: ic_put_file_model_input=0
        integer                            :: ic_put_file_model_output=0
        integer                            :: ic_get_file_model_input=0
        integer                            :: ic_get_file_model_output=0
        integer                            :: ic_put_parameter=0
        integer                            :: ic_parameter_check=0
        integer                            :: ic_get_parameter_name=0
        integer                            :: ic_get_parameter_index=0
        integer                            :: ic_process_template_files=0
        integer                            :: ic_set_working_directory=0
        integer                            :: ic_write_model_input_files=0
        integer                            :: ic_get_dimensions=0
        integer                            :: ic_get_temp_ins_status=0
        integer                            :: ic_put_observation=0
        integer                            :: ic_observation_check=0
        integer                            :: ic_get_observation_name=0
        integer                            :: ic_get_observation_index=0
        integer                            :: ic_store_instruction_set=0
        integer                            :: ic_delete_model_output_files=0
        integer                            :: ic_read_model_output_files=0
        integer                            :: ic_finalise=0

      end type MIO_STRUC
      public mio_struc

! -- Global data types are declared.

#ifdef UNIX
      integer                            :: ios=1      ! operating system
#else
      integer                            :: ios=0      ! operating system
#endif
      integer                            :: ifail
      integer                            :: ierr
      integer                            :: ip_last=0

      character (len=1)                  :: atemp1
      character (len=10)                 :: atemp10
      character (len=10)                 :: atemp12  !   MD??? error???? Should be 12
      character (len=15)                 :: atemp15
      character (len=20)                 :: atemp20
      character (len=25)                 :: atemp25
      character (len=256)                :: afile
      character (len=80)                 :: errsub
      character (len=500)                :: amessage=' '
      character (len=2000)               :: dline
      character (len=100)                :: function_name


! -- FUNCTIONS

! -- Visible functions

      public mio_initialise,                        &
             mio_initialise_model_input,            &
             mio_initialise_model_output,           &
             mio_initialise_parameters,             &
             mio_initialise_observations

      public mio_set_number_precision,              &
             mio_set_number_decpoint,               &
             mio_set_model_command,                 &
             mio_set_working_directory

      public mio_get_model_command,                 &
             mio_get_dimensions,                    &
             mio_get_temp_ins_status

      public mio_put_file_model_input,              &
             mio_put_file_model_output,             &
             mio_put_parameter,                     &
             mio_put_observation

      public mio_get_file_model_input,              &
             mio_get_file_model_output,             &
             mio_get_parameter_name,                &
             mio_get_parameter_index,               &
             mio_get_observation_name,              &
             mio_get_observation_index

      public mio_parameter_check,                   &
             mio_observation_check

      public mio_process_template_files,            &
             mio_write_model_input_files

      public mio_store_instruction_set,             &
             mio_delete_model_output_files,         &
             mio_read_model_output_files

      public mio_finalise

contains



integer function mio_initialise(estruc,mstruc)

! -- Function MIO_INITIALISE initialises the model_input_output module.

       implicit none
      
       type(err_failure_struc), intent(inout)   :: estruc     ! for error messages
       type(mio_struc), intent(inout)           :: mstruc     ! structure holding current mio dataset

       function_name='MIO_INITIALISE'
       mstruc%ic_initialise=mstruc%ic_initialise+1
       errsub='Error in function '//trim(function_name)//':'
       ifail=0

! -- Variables are initialised.

       ip_last=0
       amessage=' '

       mstruc%asize=0
       mstruc%numl=0
       mstruc%ninstr=0
       mstruc%template_status=0
       mstruc%instruction_status=0
       mstruc%precis=0
       mstruc%nopnt=0
       mstruc%mcall=0
       mstruc%workdir=' '
       mstruc%comline=' '
       mstruc%npar=0
       mstruc%nobs=0
       mstruc%numinfile=0
       mstruc%numoutfile=0

       mstruc%ic_initialise_model_input=0
       mstruc%ic_initialise_model_output=0
       mstruc%ic_initialise_parameters=0
       mstruc%ic_initialise_observations=0
       mstruc%ic_set_number_precision=0
       mstruc%ic_set_number_decpoint=0
       mstruc%ic_set_model_command=0
       mstruc%ic_get_model_command=0
       mstruc%ic_put_file_model_input=0
       mstruc%ic_put_file_model_output=0
       mstruc%ic_get_file_model_input=0
       mstruc%ic_get_file_model_output=0
       mstruc%ic_put_parameter=0
       mstruc%ic_parameter_check=0
       mstruc%ic_get_parameter_name=0
       mstruc%ic_get_parameter_index=0
       mstruc%ic_process_template_files=0
       mstruc%ic_set_working_directory=0
       mstruc%ic_write_model_input_files=0
       mstruc%ic_get_dimensions=0
       mstruc%ic_get_temp_ins_status=0
       mstruc%ic_put_observation=0
       mstruc%ic_observation_check=0
       mstruc%ic_get_observation_name=0
       mstruc%ic_get_observation_index=0
       mstruc%ic_store_instruction_set=0
       mstruc%ic_delete_model_output_files=0
       mstruc%ic_read_model_output_files=0
       mstruc%ic_finalise=0
      
       go to 9900

! -- Function exit point.

9900   continue
       if(ifail.eq.0)then    ! MD??? How ifail could be different of 0??? 
         mio_initialise=0
       else if(ifail.eq.1)then
         mio_initialise=1
         call err_reset(estruc)
         call err_add_error(estruc,amessage,function_name)
       endif
       return

end function mio_initialise



integer function mio_initialise_model_input(estruc,mstruc,numinfile)

! -- Function MIO_INITIALISE_MODEL_INPUT receives the number of template and model input files.

        implicit none
      
        type(err_failure_struc), intent(inout)   :: estruc     ! for error messages
        type(mio_struc), intent(inout)           :: mstruc     ! structure holding current mio dataset
        integer, intent(in)                      :: numinfile  ! number of model input files

        function_name='MIO_INITIALISE_MODEL_INPUT'
        mstruc%ic_initialise_model_input=mstruc%ic_initialise_model_input+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

        if(mstruc%ic_initialise.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE has been called.')
          ifail=1
          go to 9900
        endif

        if(numinfile.le.0)then
          write(amessage,10) trim(errsub)
10        format(a,' NUMINFILE variable must be supplied as positive.')
          ifail=1
          go to 9900
        endif
        allocate(mstruc%tempfile(numinfile),mstruc%modinfile(numinfile),mstruc%pardel(numinfile),stat=ierr) !MD??? Is the number of files always the same for all types?
        if(ierr.ne.0) go to 9200
        mstruc%tempfile=' '        ! an array
        mstruc%modinfile=' '       ! an array
        mstruc%pardel=' '          ! an array
        mstruc%numinfile=numinfile

        go to 9900

9200    write(amessage,9210) trim(errsub)
9210    format(a,' cannot allocate sufficient memory to store model interface work arrays.')
        ifail=1
        go to 9900
      
! -- Function exit point.

9900    continue
        if(ifail.eq.0)then
          mio_initialise_model_input=0
        else if(ifail.eq.1)then
          mio_initialise_model_input=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_initialise_model_input




integer function mio_initialise_model_output(estruc,mstruc,numoutfile)

! -- Function MIO_INITIALISE_MODEL_OUTPUT receives the number of instruction and model output files.

        implicit none
      
        type(err_failure_struc), intent(inout)   :: estruc     ! for error messages
        type(mio_struc), intent(inout)           :: mstruc     ! structure holding current mio dataset
        integer, intent(in)                      :: numoutfile ! number of model output files

        function_name='MIO_INITIALISE_MODEL_OUTPUT'
        mstruc%ic_initialise_model_output=mstruc%ic_initialise_model_output+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

        if(mstruc%ic_initialise.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE has been called.')
          ifail=1
          go to 9900
        endif

        if(numoutfile.le.0)then
          write(amessage,10) trim(errsub)
10        format(a,' NUMOUTFILE variable must be supplied as positive.')
          ifail=1
          go to 9900
        endif
        allocate(mstruc%insfile(numoutfile),mstruc%modoutfile(numoutfile),mstruc%mrkdel(numoutfile),stat=ierr)
        if(ierr.ne.0) go to 9200
        mstruc%insfile=' '         ! an array
        mstruc%modoutfile=' '      ! an array
        mstruc%mrkdel=' '          ! an array
        mstruc%numoutfile=numoutfile

        go to 9900

9200    write(amessage,9210) trim(errsub)
9210    format(a,' cannot allocate sufficient memory to store model interface work arrays.')
        ifail=1
        go to 9900
      
! -- Function exit point.

9900    continue
        if(ifail.eq.0)then
          mio_initialise_model_output=0
        else if(ifail.eq.1)then
          mio_initialise_model_output=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_initialise_model_output



integer function mio_initialise_parameters(estruc,mstruc,npar)

! -- Function MIO_INITIALISE_PARAMETERS initialises parameter storage in the model_input_output module.

        implicit none
      
        type(err_failure_struc), intent(inout)   :: estruc     ! for error messages
        type(mio_struc), intent(inout)           :: mstruc     ! structure holding current mio dataset
        integer, intent(in)                      :: npar       ! number of parameters

        function_name='MIO_INITIALISE_PARAMETERS'
        mstruc%ic_initialise_parameters=mstruc%ic_initialise_parameters+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

        if(mstruc%ic_initialise.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE has been called.')
          ifail=1
          go to 9900
        endif

        if(npar.le.0)then
          write(amessage,30) trim(errsub)
30        format(a,' NPAR variable must be supplied as positive.')
          ifail=1
          go to 9900
        endif
        mstruc%npar=npar

        allocate(mstruc%nw(npar),mstruc%pword(npar),mstruc%apar(npar),stat=ierr)
        if(ierr.ne.0) go to 9200
        mstruc%nw=0        ! an array
        mstruc%pword=' '   ! an array
        mstruc%apar=' '    ! an array

        go to 9900

9200    write(amessage,9210) trim(errsub)
9210    format(a,' cannot allocate sufficient memory to store model interface work arrays.')
        ifail=1
        go to 9900
      
! -- Function exit point.

9900    continue
        if(ifail.eq.0)then
          mio_initialise_parameters=0
        else if(ifail.eq.1)then
          mio_initialise_parameters=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_initialise_parameters



integer function mio_initialise_observations(estruc,mstruc,nobs)

! -- Function MIO_INITIALISE_OBSERVATIONS initialises observation storage in the model_input_output module.

        implicit none
      
        type(err_failure_struc), intent(inout)   :: estruc     ! for error messages
        type(mio_struc), intent(inout)           :: mstruc     ! structure holding current mio dataset
        integer, intent(in)                      :: nobs       ! number of observations
      
        function_name='MIO_INITIALISE_OBSERVATIONS'
        mstruc%ic_initialise_observations=mstruc%ic_initialise_observations+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

        if(mstruc%ic_initialise.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE has been called.')
          ifail=1
          go to 9900
        endif

        if(nobs.le.0)then
          write(amessage,30) trim(errsub)
30        format(a,' NOBS variable must be supplied as positive.')
          ifail=1
          go to 9900
        endif
        mstruc%nobs=nobs

        allocate(mstruc%obsn1(nobs),mstruc%obsn2(nobs),mstruc%iiobs(nobs),mstruc%aobs(nobs),stat=ierr)
        if(ierr.ne.0) go to 9200
        mstruc%obsn1=0
        mstruc%obsn2=0
        mstruc%iiobs=0
        mstruc%aobs=' '

        go to 9900

9200    write(amessage,9210) trim(errsub)
9210    format(a,' cannot allocate sufficient memory to store model interface work arrays.')
        ifail=1
        go to 9900
      
! -- Function exit point.

9900    continue
        if(ifail.eq.0)then
          mio_initialise_observations=0
        else if(ifail.eq.1)then
          mio_initialise_observations=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_initialise_observations



integer function mio_set_number_precision(estruc,mstruc,precision)

! -- Function MIO_SET_NUMBER_PRECISION sets the precision with which numbers are written to model
!    input files.

        implicit none
      
        type(err_failure_struc), intent(inout)   :: estruc      ! for error messages
        type(mio_struc), intent(inout)           :: mstruc      ! structure holding current mio dataset
        character (len=*), intent(in)            :: precision   ! the precision type for writing numbers

        function_name='MIO_SET_NUMBER_PRECISION'
        mstruc%ic_set_number_precision=mstruc%ic_set_number_precision+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

        if(mstruc%ic_initialise.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE has been called.')
          ifail=1
          go to 9900
        endif

        atemp12=adjustl(precision)
        call utl_casetrans(atemp12,'lo')
        if(atemp12(1:6).eq.'double')then
          mstruc%precis=1
        else if(atemp12(1:6).eq.'single')then
          mstruc%precis=0
        else
          write(amessage,30) trim(errsub)
30        format(a,' the PRECISION argument must be supplied as "single" or "double".')
          ifail=1
        endif

! -- Function exit point.

9900    continue
        if(ifail.eq.0)then
          mio_set_number_precision=0
        else if(ifail.eq.1)then
          mio_set_number_precision=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif

        return

end function mio_set_number_precision



integer function mio_set_number_decpoint(estruc,mstruc,decpoint)

! -- Function MIO_SET_NUMBER_DECPOINT sets whether a decimal point is excluded (if possible)
!    in writing a number to a model input file.

        implicit none
      
        type(err_failure_struc), intent(inout)   :: estruc      ! for error messages
        type(mio_struc), intent(inout)           :: mstruc      ! structure holding current mio dataset
        character (len=*), intent(in)            :: decpoint    ! whether to exclude a decimal point if possible

        function_name='MIO_SET_NUMBER_DECPOINT'
        mstruc%ic_set_number_decpoint=mstruc%ic_set_number_decpoint+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

        if(mstruc%ic_initialise.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE has been called.')
          ifail=1
          go to 9900
        endif

        atemp12=adjustl(decpoint)
        call utl_casetrans(atemp12,'lo')
        if(atemp12(1:7).eq.'nopoint')then
          mstruc%nopnt=1
        else if(atemp12(1:5).eq.'point')then
          mstruc%nopnt=0
        else
          write(amessage,30) trim(errsub)
30        format(a,' the DECPOINT argument must be supplied as "point" or "nopoint".')
          ifail=1
        endif

! -- Function exit point.

9900    continue
        if(ifail.eq.0)then
          mio_set_number_decpoint=0
        else if(ifail.eq.1)then
          mio_set_number_decpoint=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif

        return

end function mio_set_number_decpoint



integer function mio_set_model_command(estruc,mstruc,command)

! -- Function MIO_SET_MODEL_COMMAND sets the command to run the model.

        implicit none
      
        type(err_failure_struc), intent(inout)   :: estruc      ! for error messages
        type(mio_struc), intent(inout)           :: mstruc      ! structure holding current mio dataset
        character (len=*), intent(in)            :: command     ! command to run model

        function_name='MIO_SET_MODEL_COMMAND'
        mstruc%ic_set_model_command=mstruc%ic_set_model_command+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

        if(mstruc%ic_initialise.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE has been called.')
          ifail=1
          go to 9900
        endif

        mstruc%comline=adjustl(command)
        if(ios.eq.0) call utl_casetrans(mstruc%comline,'lo')

9900    continue
        if(ifail.eq.0)then
          mio_set_model_command=0
        else if(ifail.eq.1)then
          mio_set_model_command=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif

        return

end function mio_set_model_command



integer function mio_get_model_command(estruc,mstruc,command)

! -- Function MIO_GET_MODEL_COMMAND retreives the command to run the model.

        implicit none
      
        type(err_failure_struc), intent(inout)   :: estruc      ! for error messages
        type(mio_struc), intent(inout)           :: mstruc      ! structure holding current mio dataset
        character (len=*), intent(out)           :: command     ! command to run model

        integer                                  :: n

        function_name='MIO_GET_MODEL_COMMAND'
        mstruc%ic_get_model_command=mstruc%ic_get_model_command+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

        if(mstruc%ic_initialise.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE has been called.')
          ifail=1
          go to 9900
        endif

        n=len_trim(mstruc%comline)
        if(n.gt.len(command)) n=len(command)
        command=mstruc%comline(1:n)

9900    continue
        if(ifail.eq.0)then
          mio_get_model_command=0
        else if(ifail.eq.1)then
          mio_get_model_command=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif

        return

end function mio_get_model_command




integer function mio_put_file_model_input(estruc,mstruc,inum,filename1,filename2)

! -- Function MIO_PUT_FILE_MODEL_INPUT supplies the name of a template and model input file.

        implicit none
       
        type(err_failure_struc), intent(inout)     :: estruc      ! for error messages
        type(mio_struc), intent(inout)             :: mstruc      ! structure holding current mio dataset
        integer, intent(in)                        :: inum        ! file number
        character (len=*), intent(in)              :: filename1   ! name of template file
        character (len=*), intent(in)              :: filename2   ! name of model input file
       
        function_name='MIO_PUT_FILE_MODEL_INPUT'
        mstruc%ic_put_file_model_input=mstruc%ic_put_file_model_input+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

        if(mstruc%ic_initialise_model_input.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE_MODEL_INPUT has been called.')
          ifail=1
          go to 9900
        endif

        if((inum.lt.1).or.(inum.gt.mstruc%numinfile))then
          write(amessage,10) trim(errsub)
10        format(a,' INUM argument out of range.')
          ifail=1
          go to 9900
        endif

        afile=adjustl(filename1)
        if(ios.eq.0) call utl_casetrans(afile,'lo')
        mstruc%tempfile(inum)=afile
        afile=adjustl(filename2)
        if(ios.eq.0) call utl_casetrans(afile,'lo')
        mstruc%modinfile(inum)=afile

        go to 9900

! -- Function exit point.

9900    continue
        if(ifail.eq.0)then
          mio_put_file_model_input=0
        else if(ifail.eq.1)then
          mio_put_file_model_input=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_put_file_model_input




integer function mio_put_file_model_output(estruc,mstruc,inum,filename1,filename2)

! -- Function MIO_PUT_FILE_MODEL_OUTPUT supplies the name of an instruction and model output file.

        implicit none
       
        type(err_failure_struc), intent(inout)     :: estruc      ! for error messages
        type(mio_struc), intent(inout)             :: mstruc      ! structure holding current mio dataset
        integer, intent(in)                        :: inum        ! file number
        character (len=*), intent(in)              :: filename1   ! name of instruction file
        character (len=*), intent(in)              :: filename2   ! name of model output file

        integer                                    :: icount,i,j,k
       
        function_name='MIO_PUT_FILE_MODEL_OUTPUT'
        mstruc%ic_put_file_model_output=mstruc%ic_put_file_model_output+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

        if(mstruc%ic_initialise_model_output.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE_MODEL_OUTPUT has been called.')
          ifail=1
          go to 9900
        endif

        if((inum.lt.1).or.(inum.gt.mstruc%numoutfile))then
          write(amessage,10) trim(errsub)
10        format(a,' INUM argument out of range.')
          ifail=1
          go to 9900
        endif

        afile=adjustl(filename1)
        if(ios.eq.0) call utl_casetrans(afile,'lo')
        mstruc%insfile(inum)=afile
        afile=adjustl(filename2)
        if(ios.eq.0) call utl_casetrans(afile,'lo')
        mstruc%modoutfile(inum)=afile
        if(mstruc%asize.ne.0)then
          icount=0
          do i=1,mstruc%asize
            if(mstruc%a(i).eq.achar(2))then
              icount=icount+1
              if(icount.eq.inum)then
                k=0
                do j=i+2,i+2+len(mstruc%modoutfile(inum))-1
                  k=k+1
                  if(k.gt.len(mstruc%modoutfile(inum)))exit
                  mstruc%a(j)=mstruc%modoutfile(inum)(k:k)
                enddo
                go to 8000
              endif
            endif
          enddo
        endif

8000    continue

        go to 9900

! -- Function exit point.

9900    continue
        if(ifail.eq.0)then
          mio_put_file_model_output=0
        else if(ifail.eq.1)then
          mio_put_file_model_output=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_put_file_model_output



integer function mio_get_file_model_input(estruc,mstruc,inum,filename1,filename2)

! -- Function MIO_GET_FILE_MODEL_INPUT retreives the name of a template and model input file.

        implicit none
       
        type(err_failure_struc), intent(inout)     :: estruc      ! for error messages
        type(mio_struc), intent(inout)             :: mstruc      ! structure holding current mio dataset
        integer, intent(in)                        :: inum        ! file number
        character (len=*), intent(out)             :: filename1   ! name of template file
        character (len=*), intent(out)             :: filename2   ! name of model input file


        function_name='MIO_GET_FILE_MODEL_INPUT'
        mstruc%ic_get_file_model_input=mstruc%ic_get_file_model_input+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

        if(mstruc%ic_initialise_model_input.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE_MODEL_INPUT has been called.')
          ifail=1
          go to 9900
        endif

        if((inum.lt.1).or.(inum.gt.mstruc%numinfile))then
          write(amessage,10) trim(errsub)
10        format(a,' INUM argument out of range.')
          ifail=1
          go to 9900
        endif
        filename1=mstruc%tempfile(inum)
        filename2=mstruc%modinfile(inum)

9900    continue
        if(ifail.eq.0)then
          mio_get_file_model_input=0
        else if(ifail.eq.1)then
          mio_get_file_model_input=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_get_file_model_input



integer function mio_get_file_model_output(estruc,mstruc,inum,filename1,filename2)

! -- Function MIO_GET_FILE OUTPUT retrieves the name of an instruction and model output file.

        implicit none
       
        type(err_failure_struc), intent(inout)     :: estruc      ! for error messages
        type(mio_struc), intent(inout)             :: mstruc      ! structure holding current mio dataset
        integer, intent(in)                        :: inum        ! file number
        character (len=*), intent(out)             :: filename1   ! name of instruction
        character (len=*), intent(out)             :: filename2   ! name of model output file

        function_name='MIO_GET_FILE_MODEL_OUTPUT'
        mstruc%ic_get_file_model_output=mstruc%ic_get_file_model_output+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

        if(mstruc%ic_initialise_model_output.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE_MODEL_OUTPUT has been called.')
          ifail=1
          go to 9900
        endif

        if((inum.lt.1).or.(inum.gt.mstruc%numoutfile))then
          write(amessage,10) trim(errsub)
10        format(a,' INUM argument out of range.')
          ifail=1
          go to 9900
        endif
        filename1=mstruc%insfile(inum)
        filename2=mstruc%modoutfile(inum)

9900    continue
        if(ifail.eq.0)then
          mio_get_file_model_output=0
        else if(ifail.eq.1)then
          mio_get_file_model_output=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_get_file_model_output




integer function mio_put_parameter(estruc,mstruc,inum,apar)

! -- Function MIO_PUT_PARAMETER supplies the name of a parameter.

        implicit none
       
        type(err_failure_struc), intent(inout)     :: estruc      ! for error messages
        type(mio_struc), intent(inout)             :: mstruc      ! structure holding current mio dataset
        integer, intent(in)                        :: inum        ! parameter number
        character (len=*), intent(in)              :: apar        ! name of parameter

        integer                                    :: n1
       
        function_name='MIO_PUT_PARAMETER'
        mstruc%ic_put_parameter=mstruc%ic_put_parameter+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

! -- Context checking is undertaken.

        if(mstruc%ic_initialise_parameters.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE_PARAMETERS has been called.')
          ifail=1
          go to 9900
        endif

! -- Other checking is undertaken.

        if((inum.lt.1).or.(inum.gt.mstruc%npar)) then
          write(amessage,10) trim(errsub)
10        format(a,' supplied parameter index out of range.')
          ifail=1
          go to 9900
        endif
        atemp15=adjustl(apar)
        n1=len_trim(atemp15)
        if(n1.gt.12)then
          write(amessage,20) trim(errsub),trim(atemp15)
20        format(a,' supplied parameter name "',a,'" greater than 12 characters in length.')
          ifail=1
          go to 9900
        endif
        call utl_casetrans(atemp15,'lo')
        mstruc%apar(inum)=trim(atemp15)

9900    continue
        if(ifail.eq.0)then
          mio_put_parameter=0
        else if(ifail.eq.1)then
          mio_put_parameter=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_put_parameter



integer function mio_parameter_check(estruc,mstruc)

! -- Function MIO_PARAMETER_CHECK checks the integrity of a supplied parameter set.
!    Note that it returns on encountering the first error.

        implicit none
       
        type(err_failure_struc), intent(inout)     :: estruc      ! for error messages
        type(mio_struc), intent(inout)             :: mstruc      ! structure holding current mio dataset

        integer                                    :: ip,np,n2,i
        character (len=12)                         :: aword
        character (len=12), allocatable            :: awork(:)

        function_name='MIO_PARAMETER_CHECK'
        mstruc%ic_parameter_check=mstruc%ic_parameter_check+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

! -- Context checking is undertaken.

        if(mstruc%ic_initialise_parameters.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE_PARAMETERS has been called.')
          ifail=1
          go to 9900
        endif

! -- Individual parameter names are checked.

        errsub = 'Error in model parameter dataset:'
        np=mstruc%npar
        do ip=1,np
          atemp12=mstruc%apar(ip)
          if(atemp12.eq.' ')then
            call utl_int2char(ip,atemp10)
            write(amessage,20) trim(errsub),trim(atemp10)
20          format(a,' no value has been supplied for parameter number ',a,'.')
            go to 9800
          endif
          n2=index(trim(atemp12),' ')
          if(n2.ne.0)then
            write(amessage,30) trim(errsub),trim(atemp12)
30          format(a,' parameter name "',a,'" contains a space.')
            go to 9800
          endif
          if(index(atemp12,char(9)).ne.0)then
            write(amessage,40) trim(errsub),trim(atemp12)
40          format(a,' parameter name "',a,'" contains a tab.')
            go to 9800
          endif
        enddo

! -- Duplicate parameter names are checked for.

        allocate(awork(np),stat=ierr)
        if(ierr.ne.0)then
          write(amessage,50)
50        format('Error in function MIO_PARAMETER_CHECK: ',  &
          ' cannot allocate sufficient memory to continue execution.')
          go to 9800
        endif
        i=utl_duplicate_check(np,mstruc%apar,awork,aword)
        deallocate(awork,stat=ierr)
        if(i.ne.0)then
          write(amessage,60) trim(errsub),trim(aword)
60        format(a,' parameter name "',a,'" duplicated in supplied list of parameters.')
          go to 9800
        endif

! -- The presence of all template and model input files is checked for.

        do i=1,mstruc%numinfile
          if(mstruc%tempfile(i).eq.' ')then
            call utl_int2char(i,atemp10)
            write(amessage,80) trim(errsub),trim(atemp10)
80          format(a,' a name has not been provided for template file number ',a,'.')
            go to 9800
          endif
          if(mstruc%modinfile(i).eq.' ')then
            call utl_int2char(i,atemp10)
            write(amessage,90) trim(errsub),trim(atemp10)
90          format(a,' a name has not been provided for model input file number ',a,'.')
            go to 9800
          endif
        enddo

        go to 9900

9800    ifail=1
9900    continue
        if(ifail.eq.0)then
          mio_parameter_check=0
        else if(ifail.eq.1)then
          mio_parameter_check=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_parameter_check



integer function mio_get_parameter_name(estruc,mstruc,inum,apar)

! -- Function MIO_GET_PARAMETER_NAME retreives the name of a parameter, given its index.

        implicit none

        type(err_failure_struc), intent(inout)     :: estruc
        type(mio_struc), intent(inout)             :: mstruc
        integer, intent(in)                        :: inum
        character (len=*), intent(out)             :: apar

        function_name='MIO_GET_PARAMETER_NAME'
        mstruc%ic_get_parameter_name=mstruc%ic_get_parameter_name+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

        if(mstruc%ic_initialise_parameters.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE_PARAMETERS has been called.')
          ifail=1
          go to 9900
        endif

        if((inum.lt.1).or.(inum.gt.mstruc%npar)) then
          write(amessage,10) trim(errsub)
10        format(a,' supplied parameter index out of range.')
          ifail=1
          go to 9900
        endif

        apar=mstruc%apar(inum)

9900    continue
        if(ifail.eq.0)then
          mio_get_parameter_name=0
        else if(ifail.eq.1)then
          mio_get_parameter_name=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_get_parameter_name



integer function mio_get_parameter_index(estruc,mstruc,apar,inum)

! -- Function MIO_GET_PARAMETER_INDEX retreives the index of a parameter, given its name.

        implicit none

        type(err_failure_struc), intent(inout)     :: estruc
        type(mio_struc), intent(inout)             :: mstruc
        character (len=*), intent(in)              :: apar
        integer, intent(out)                       :: inum

        integer                                    :: ip,np,jfail

        function_name='MIO_GET_PARAMETER_INDEX'
        mstruc%ic_get_parameter_index=mstruc%ic_get_parameter_index+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0
        inum=0

        if(mstruc%ic_initialise_parameters.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE_PARAMETERS has been called.')
          ifail=1
          go to 9900
        endif

        atemp15=adjustl(apar)
        ip=len_trim(atemp15)
        if(ip.gt.12)then
          write(amessage,20) trim(errsub),trim(atemp15)
20        format(a,' supplied parameter name "',a,'" is greater than 12 characters in length.')
          ifail=1
          go to 9900
        endif

        ip=ip_last
        atemp12=adjustl(apar)
        call utl_casetrans(atemp12,'lo')
        np=mstruc%npar
        call mio_which1(jfail,np,ip,mstruc%apar,atemp12)
        if(jfail.ne.0)then
          write(amessage,10) trim(errsub),trim(atemp12)
10        format(a,' parameter name "',a,'" not found in parameter list.')
          ifail=1
        else
          ip_last=ip
          inum=ip
        endif

9900    continue
        if(ifail.eq.0)then
          mio_get_parameter_index=0
        else if(ifail.eq.1)then
          mio_get_parameter_index=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_get_parameter_index




integer function mio_process_template_files(estruc,mstruc)

! -- Function MIO_PROCESS_TEMPLATE_FILES does rudmentary checking of template files.
!    However its main role is to find the smallest character width to which each
!    parameter will be written.

        implicit none
        
        type(err_failure_struc), intent(inout)     :: estruc
        type(mio_struc), intent(inout)             :: mstruc

        integer                                    :: ipar,i,iline,nblc,j2,j1,jfail,nnw,iunit
        character (len=10)                         :: aline
        character (len=12)                         :: tpar
        
        function_name='MIO_PROCESS_TEMPLATE_FILES'
        mstruc%ic_process_template_files=mstruc%ic_process_template_files+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

! -- Context checking is undertaken.

        if(mstruc%ic_parameter_check.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_PARAMETER_CHECK has been called.')
          ifail=1
          go to 9900
        endif
        if(mstruc%ic_process_template_files.eq.2)then
          mstruc%ic_process_template_files=1
          write(amessage,10) trim(errsub)
10        format(a,' this function has been called already and should not be called again ',   &
          'without re-initialisation of MODEL_INPUT_OUTPUT module.')
          ifail=1
          go to 9900
        endif

! -- Processing now begins.

        ipar=1
        do i=1,mstruc%npar
          mstruc%nw(i)=1000
        enddo

        do i=1,mstruc%numinfile
          call utl_addquote(mstruc%tempfile(i),afile)
          iunit=utl_nextunit()
          open(unit=iunit,file=mstruc%tempfile(i),status='old',iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,410) trim(afile)
410         format('Cannot open template file ',a,'.')
            ifail=1
            go to 9900
          endif
          read(iunit,'(a)',err=9000,end=9200) dline
          call utl_casetrans(dline(1:3),'lo')
          if((dline(1:3).ne.'ptf').and.(dline(1:3).ne.'jtf'))go to 9200
          mstruc%pardel(i)=dline(5:5)
          if(mstruc%pardel(i).eq.' ') go to 9200
          iline=1
520       iline=iline+1
          read(iunit,'(a)',err=9000,end=680) dline
          nblc=len_trim(dline)
          j2=0
550       if(j2.ge.nblc) go to 520
          j1=index(dline(j2+1:nblc),mstruc%pardel(i))
          if(j1.eq.0) go to 520
          j1=j1+j2
          j2=index(dline(j1+1:nblc),mstruc%pardel(i))
          if(j2.eq.0)then
            call utl_int2char(iline,aline)
            write(amessage,555) trim(aline),trim(afile)
555         format('Unbalanced parameter delimiters at line ',a,' of template file ',a,'.')
            ifail=1
            go to 9900
          endif
          j2=j2+j1
          call mio_parnam(jfail,j1,j2,tpar)
          if(jfail.eq.1)then
            call utl_int2char(iline,aline)
            write(amessage,556) trim(aline),trim(afile)
556         format('Parameter space less than three characters wide at line ',a,  &
            ' of file ',a,'.')
            ifail=1
            go to 9900
          else if (jfail.eq.2)then
            call utl_int2char(iline,aline)
            write(amessage,557) trim(aline),trim(afile)
557         format('Blank parameter space at line ',a,' of file ',a,'.')
            ifail=1
            go to 9900
          endif
          call mio_which1(jfail,mstruc%npar,ipar,mstruc%apar,tpar)
          if(jfail.ne.0)then
            call utl_int2char(iline,aline)
            write(amessage,558) trim(tpar),trim(aline),trim(afile)
558         format('Parameter "',a,'" cited on line ',a,' of template file ',a,   &
            ' has not been cited in parameter list.')
            ifail=1
            go to 9900
          endif
          nnw=j2-j1+1
          if(nnw.lt.mstruc%nw(ipar)) mstruc%nw(ipar)=nnw
          go to 550
680       close(unit=iunit)
        enddo

        mstruc%template_status=1

        go to 9900

9000    write(amessage,9010) trim(afile)
9010    format('Unable to read template file ',a,'.')
        ifail=1
        go to 9900
9200    write(amessage,9210) trim(afile)
9210    format('"ptf" or "jtf" header, followed by space, followed by parameter delimiter ',  &
        'expected on first line of template file ',a,'.')
        ifail=1
        go to 9900

9900    continue
        if(ifail.eq.0)then
          mio_process_template_files=0
        else if(ifail.eq.1)then
          mio_process_template_files=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_process_template_files




integer function mio_set_working_directory(estruc,mstruc,workingdir)

! -- Function MIO_SET_WORKING_DIRECTORY supplies a subdirectory whose name will be affixed to
!    that of all input and output files.

! -- Note that a forward or back slash is appended if appropriate.

        implicit none

        type(err_failure_struc), intent(inout)   :: estruc
        type(mio_struc), intent(inout)           :: mstruc
        character (len=*), intent(in)            :: workingdir

        integer                                  :: n

        function_name='MIO_SET_WORKING_DIRECTORY'
        mstruc%ic_set_working_directory=mstruc%ic_set_working_directory+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

        if(mstruc%ic_initialise.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE has been called.')
          ifail=1
          go to 9900
        endif

        afile=workingdir
        if(ios.eq.0) call utl_casetrans(afile,'lo')

        n=len_trim(afile)
        if(n.ne.0)then
          if(ios.eq.1)then
            atemp1='\'
          else
            atemp1='/'
          endif
          if(afile(n:n).eq.atemp1)then
            write(amessage,10) trim(errsub),trim(afile)
10          format(a,' illegal final character in supplied working directory name "',a,'".')
            ifail=1
            go to 9900
          endif

          if(ios.eq.1)then
            atemp1='/'
          else
            atemp1='\'
          endif
          if(afile(n:n).ne.atemp1)afile(n+1:n+1)=atemp1
        endif
        mstruc%workdir=trim(afile)

9900    continue
        if(ifail.eq.0)then
          mio_set_working_directory=0
        else if(ifail.eq.1)then
          mio_set_working_directory=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_set_working_directory




integer function mio_write_model_input_files(estruc,mstruc,pval)

! -- Function MIO_WRITE_MODEL_INPUT_FILES writes model input files based on a set of model
!    input template files and a set of current parameter values.

        implicit none
        
        type(err_failure_struc), intent(inout)               :: estruc
        type(mio_struc), intent(inout)                       :: mstruc
        double precision, intent(inout), dimension(:)        :: pval

        integer                                              :: ipar,ipp,jfail,ifile,iunit,iunit1,iline, &
                                                                lc,j1,j2,j,idir,n
        double precision                                     :: tval
        character (len=1)                                    :: aa
        character (len=12)                                   :: tpar
        character (len=256)                                  :: mifile

        function_name='MIO_WRITE_MODEL_INPUT_FILES'
        mstruc%ic_write_model_input_files=mstruc%ic_write_model_input_files+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

! -- Preliminary checking is undertaken.

        if(mstruc%ic_process_template_files.eq.0)then
          write(amessage,3) trim(errsub)
3         format(a,' this function must not be called unless function MIO_PROCESS_TEMPLATE_FILES ',  &
          'has been called previously.')
          ifail=1
          go to 9900
        endif

        if(size(pval).lt.mstruc%npar)then
          write(amessage,5) trim(errsub)
5         format(a,' size of PVAL array supplied to function MIO_WRITE_MODEL_INPUT_FILES is insufficient.')
          ifail=1
          go to 9900
        endif

! -- The writing is now actually done.

        errsub='Error writing parameter values to model input files:'
        idir=0
        if(mstruc%workdir.ne.' ') idir=1

! -- Each of the parameter words is filled.

        ipar=1
        do 100 ipp=1,mstruc%npar
          call mio_wrtsig(jfail,pval(ipp),mstruc%pword(ipp),mstruc%nw(ipp),mstruc%precis,tval,mstruc%nopnt)
          if(jfail.lt.0)then
            write(amessage,10) trim(errsub),trim(mstruc%apar(ipp))
10          format(a,' internal error condition has arisen while attempting to write ', &
            'current value of parameter "',a,'" to model input file.')
            ifail=1
            go to 9900
          else if (jfail.eq.1)then
            write(amessage,11) trim(errsub),trim(mstruc%apar(ipp))
11          format(a,' exponent of parameter "',a,'" is too large or too small for ', &
            'single precision protocol.')
            ifail=1
            go to 9900
          else if (jfail.eq.2)then
            write(amessage,12) trim(errsub),trim(mstruc%apar(ipp))
12          format(a,' exponent of parameter "',a,'" is too large or too small for ', &
            'double precision protocol.')
            ifail=1
            go to 9900
          else if (jfail.eq.3)then
            write(amessage,13) trim(errsub),trim(mstruc%apar(ipp))
13          format(a,' field width of parameter "',a,'" on at least one template file ', &
            'is too small to represent current parameter value. The number is too large ', &
            'to fit, or too small to be represented with any precision.')
            ifail=1
            go to 9900
          endif
          pval(ipp)=tval
100     continue

! -- Next the substitutions in the template files are made.

        do 500 ifile=1,mstruc%numinfile
          call utl_addquote(mstruc%tempfile(ifile),afile)
          iunit=utl_nextunit()
          open(unit=iunit,file=mstruc%tempfile(ifile),status='old',iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,110) trim(errsub),trim(afile)
110         format(a,' cannot open template file ',a,'.')
            ifail=1
            go to 9900
          endif
          iunit1=utl_nextunit()
          if(idir.eq.0)then
            mifile=mstruc%modinfile(ifile)
          else
            mifile=trim(mstruc%workdir)//trim(mstruc%modinfile(ifile))
          endif
          open(unit=iunit1,file=mifile,iostat=ierr)
          if(ierr.ne.0)then
            call utl_addquote(mifile,afile)
            write(amessage,115) trim(errsub),trim(afile)
115         format(a,' cannot open model input file ',a,' to write updated parameter ',  &
            'values prior to running model.')
            ifail=1
            go to 9900
          endif
          read(iunit,*,err=9000)
          iline=1
120       iline=iline+1
          read(iunit,22,end=400,err=9000) dline
22        format(a)
          lc=len_trim(dline)
          j2=0
150       if(j2.ge.lc) go to 300
          j1=index(dline(j2+1:lc),mstruc%pardel(ifile))
          if(j1.eq.0) go to 300
          j1=j1+j2
          j2=index(dline(j1+1:lc),mstruc%pardel(ifile))
          j2=j2+j1
          call mio_parnam(jfail,j1,j2,tpar)
          call mio_which1(jfail,mstruc%npar,ipar,mstruc%apar,tpar)
!       The following works when space bigger than pword(:nblnk(pword))
!       dline(j1:j2)=pword(ipar)(:nblnk(pword(ipar)))
          do 160 j=j1,j2
            dline(j:j)=' '
160       continue
          j=len_trim(mstruc%pword(ipar))
          dline(j2-j+1:j2)=mstruc%pword(ipar)(1:j)
          go to 150
300       write(iunit1,22,err=320) trim(dline)
          go to 120
320       call utl_addquote(mifile,afile)
          write(amessage,321) trim(errsub),trim(afile)
321       format(a,' cannot write to model input file ',a,'.')
          ifail=1
          go to 9900
400       close(unit=iunit)
          close(unit=iunit1,iostat=ierr)
          if(ierr.ne.0)then
            call utl_addquote(mifile,afile)
            write(amessage,321) trim(errsub),trim(afile)
            ifail=1
            go to 9900
          endif
500     continue
        go to 9900

9000    write(amessage,9010) trim(afile)
9010    format('Unable to read template file ',a,'.')
        ifail=1
        go to 9900

! -- Function exit point.

9900    continue
        if(ifail.eq.0)then
          mio_write_model_input_files=0
        else if(ifail.eq.1)then
          mio_write_model_input_files=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_write_model_input_files




integer function mio_get_dimensions(estruc,mstruc,numinfile,numoutfile,npar,nobs)

! -- Function MIO_GET_DIMENSIONS allows a calling program to determine the number of
!    dimensions of various aspects of the process of writing/reading model input/output files.

        implicit none
        
        type(err_failure_struc), intent(inout)               :: estruc
        type(mio_struc), intent(inout)                       :: mstruc
        integer, intent(out)                                 :: numinfile
        integer, intent(out)                                 :: numoutfile
        integer, intent(out)                                 :: npar
        integer, intent(out)                                 :: nobs

        function_name='MIO_GET_DIMENSIONS'
        mstruc%ic_get_dimensions=mstruc%ic_get_dimensions+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

! -- Context checking is undertaken.

        if(mstruc%ic_initialise.eq.0)then
          write(amessage,10) trim(errsub)
10        format(a,' this function must not be called before function MIO_INITIALISE has been called.')
          ifail=1
          go to 9900
        endif

! -- Values are assigned to output variables.

        numinfile=mstruc%numinfile
        numoutfile=mstruc%numoutfile
        npar=mstruc%npar
        nobs=mstruc%nobs

! -- The function exit point.

9900    continue
        if(ifail.eq.0)then
          mio_get_dimensions=0
        else if(ifail.eq.1)then
          mio_get_dimensions=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_get_dimensions



integer function mio_get_temp_ins_status(estruc,mstruc,template_status,instruction_status)

! -- Function MIO_GET_TEMP_INS_STATUS retrieves the status of template and instruction file processing.

        implicit none

        type(err_failure_struc), intent(inout)     :: estruc              ! for error messages
        type(mio_struc), intent(inout)             :: mstruc              ! structure holding current mio dataset
        integer, intent(out)                       :: template_status     ! status of template file processing
        integer, intent(out)                       :: instruction_status  ! status of instruction file processing

        function_name='MIO_GET_TEMP_INS_STATUS'
        mstruc%ic_get_temp_ins_status=mstruc%ic_get_temp_ins_status+1
        template_status=mstruc%template_status
        instruction_status=mstruc%instruction_status

        mio_get_temp_ins_status=0
        return

end function mio_get_temp_ins_status




integer function mio_put_observation(estruc,mstruc,inum,aobs)

! -- Function MIO_PUT_OBSERVATION supplies the name of an observation.

        implicit none
       
        type(err_failure_struc), intent(inout)     :: estruc      ! for error messages
        type(mio_struc), intent(inout)             :: mstruc      ! structure holding current mio dataset
        integer, intent(in)                        :: inum        ! parameter number
        character (len=*), intent(in)              :: aobs        ! name of observation

        integer                                    :: n1
       
        function_name='MIO_PUT_OBSERVATION'
        mstruc%ic_put_observation=mstruc%ic_put_observation+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

! -- Context checking is undertaken.

        if(mstruc%ic_initialise_observations.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE_OBSERVATIONS has been called.')
          ifail=1
          go to 9900
        endif

! -- Other checking is undertaken.

        if((inum.lt.1).or.(inum.gt.mstruc%nobs)) then
          write(amessage,10) trim(errsub)
10        format(a,' supplied observation index out of range.')
          ifail=1
          go to 9900
        endif
        atemp25=adjustl(aobs)
        n1=len_trim(atemp25)
        if(n1.gt.20)then
          write(amessage,20) trim(errsub),trim(atemp25)
20        format(a,' supplied observation name "',a,'" greater than 20 characters in length.')
          ifail=1
          go to 9900
        endif
        call utl_casetrans(atemp25,'lo')
        mstruc%aobs(inum)=trim(atemp25)

9900    continue
        if(ifail.eq.0)then
          mio_put_observation=0
        else if(ifail.eq.1)then
          mio_put_observation=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_put_observation




integer function mio_observation_check(estruc,mstruc)

! -- Function MIO_OBSERVATION_CHECK checks the integrity of a supplied observation set.
!    Note that it returns immediately on encountering the first error.

        implicit none
       
        type(err_failure_struc), intent(inout)     :: estruc      ! for error messages
        type(mio_struc), intent(inout)             :: mstruc      ! structure holding current mio dataset

        integer                                    :: ip,np,n2,i
        character (len=20)                         :: aword
        character (len=20), allocatable            :: awork(:)

        function_name='MIO_OBSERVATION_CHECK'
        mstruc%ic_observation_check=mstruc%ic_observation_check+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

! -- Context checking is undertaken.

        if(mstruc%ic_initialise_observations.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE_OBSERVATIONS has been called.')
          ifail=1
          go to 9900
        endif

! -- Individual names are checked for.

        errsub = 'Error in model observation dataset:'
        np=mstruc%nobs
        do ip=1,np
          atemp20=mstruc%aobs(ip)
          if(atemp20.eq.' ')then
            call utl_int2char(ip,atemp10)
            write(amessage,20) trim(errsub),trim(atemp10)
20          format(a,' no value has been supplied for observation number ',a,'.')
            go to 9800
          endif
          n2=index(trim(atemp20),' ')
          if(n2.ne.0)then
            write(amessage,30) trim(errsub),trim(atemp20)
30          format(a,' observation name "',a,'" contains a space.')
            go to 9800
          endif
          if(index(atemp20,char(9)).ne.0)then
            write(amessage,40) trim(errsub),trim(atemp20)
40          format(a,' observation name "',a,'" contains a tab.')
            go to 9800
          endif
        enddo

! -- Duplicates are checked for.

        allocate(awork(np),stat=ierr)
        if(ierr.ne.0)then
          write(amessage,50)
50        format('Error in function MIO_OBSERVATION_CHECK: ',  &
          ' cannot allocate sufficient memory to continue execution.')
          go to 9800
        endif
        i=utl_duplicate_check(np,mstruc%aobs,awork,aword)
        deallocate(awork,stat=ierr)
        if(i.ne.0)then
          write(amessage,60) trim(errsub),trim(aword)
60        format(a,' observation name "',a,'" duplicated in supplied list of observations.')
          go to 9800
        endif

! -- The presence of all instruction and model output files is checked for.

        do i=1,mstruc%numoutfile
          if(mstruc%insfile(i).eq.' ')then
            call utl_int2char(i,atemp10)
            write(amessage,10) trim(errsub),trim(atemp10)
10          format(a,' a name has not been provided for instruction file number ',a,'.')
            go to 9800
          endif
          if(mstruc%modoutfile(i).eq.' ')then
            call utl_int2char(i,atemp10)
            write(amessage,11) trim(errsub),trim(atemp10)
11          format(a,' a name has not been provided for model output file number ',a,'.')
            go to 9800
          endif
        enddo

        go to 9900

9800    ifail=1
9900    continue
        if(ifail.eq.0)then
          mio_observation_check=0
        else if(ifail.eq.1)then
          mio_observation_check=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_observation_check



integer function mio_get_observation_name(estruc,mstruc,inum,aobs)

! -- Function MIO_GET_OBSERVATION_NAME retreives the name of an observation, given its index.

        implicit none

        type(err_failure_struc), intent(inout)     :: estruc
        type(mio_struc), intent(inout)             :: mstruc
        integer, intent(in)                        :: inum
        character (len=*), intent(out)             :: aobs

        function_name='MIO_GET_OBSERVATION_NAME'
        mstruc%ic_get_observation_name=mstruc%ic_get_observation_name+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

        if(mstruc%ic_initialise_observations.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE_OBSERVATIONS has been called.')
          ifail=1
          go to 9900
        endif

        if((inum.lt.1).or.(inum.gt.mstruc%nobs)) then
          write(amessage,10) trim(errsub)
10        format(a,' supplied observation index out of range.')
          ifail=1
          go to 9900
        endif

        aobs=mstruc%aobs(inum)

9900    continue
        if(ifail.eq.0)then
          mio_get_observation_name=0
        else if(ifail.eq.1)then
          mio_get_observation_name=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_get_observation_name



integer function mio_get_observation_index(estruc,mstruc,aobs,inum)

! -- Function MIO_GET_OBSERVATION_INDEX retreives the index of an observation, given its name.

        implicit none

        type(err_failure_struc), intent(inout)     :: estruc
        type(mio_struc), intent(inout)             :: mstruc
        character (len=*), intent(in)              :: aobs
        integer, intent(out)                       :: inum

        integer                                    :: ip,np,jfail

        function_name='MIO_GET_OBSERVATION_INDEX'
        mstruc%ic_get_observation_index=mstruc%ic_get_observation_index+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0
        inum=0

        if(mstruc%ic_initialise_observations.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_INITIALISE_OBSERVATIONS has been called.')
          ifail=1
          go to 9900
        endif

        atemp25=adjustl(aobs)
        ip=len_trim(atemp25)
        if(ip.gt.20)then
          write(amessage,20) trim(errsub),trim(atemp25)
20        format(a,' supplied observation name "',a,'" greater than 20 characters in length.')
          ifail=1
          go to 9900
        endif

        ip=ip_last
        atemp20=adjustl(aobs)
        call utl_casetrans(atemp20,'lo')
        np=mstruc%nobs
        call mio_which1(jfail,np,ip,mstruc%aobs,atemp20)
        if(jfail.ne.0)then
          write(amessage,10) trim(errsub),trim(atemp20)
10        format(a,' observation name "',a,'" not found in observation list.')
          ifail=1
        else
          ip_last=ip
          inum=ip
        endif

9900    continue
        if(ifail.eq.0)then
          mio_get_observation_index=0
        else if(ifail.eq.1)then
          mio_get_observation_index=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_get_observation_index



integer function mio_store_instruction_set(estruc,mstruc)

! -- Function MIO_STORE_INSTRUCTION_SET reads all instruction files, storing the
!    instructions contained therein for more efficient later access.

        implicit none

        type(err_failure_struc), intent(inout)          :: estruc
        type(mio_struc), intent(inout)                  :: mstruc

        integer                                         :: i,nblbmx,iunit,nblc,j,ins,isum

        function_name='MIO_STORE_INSTRUCTION_SET'
        mstruc%ic_store_instruction_set=mstruc%ic_store_instruction_set+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

! -- Context checking is undertaken.

        if(mstruc%ic_observation_check.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_OBSERVATION_CHECK has been called.')
          ifail=1
          go to 9900
        endif
        if(mstruc%ic_store_instruction_set.eq.2)then
          mstruc%ic_store_instruction_set=1
          write(amessage,10) trim(errsub)
10        format(a,' this function has been called already and should not be called again ',   &
          'without re-initialisation of MODEL_INPUT_OUTPUT module.')
          ifail=1
          go to 9900
        endif

! -- Processing is now undertaken.

        errsub='Error in instructions to read model output file(s):'
        nblbmx=0
        mstruc%asize=0
        mstruc%numl=0
        mstruc%ninstr=0
        do i=1,mstruc%numoutfile
          call utl_addquote(mstruc%insfile(i),afile)
          iunit=utl_nextunit()
          open(unit=iunit,file=mstruc%insfile(i),status='old',iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,20) trim(errsub),trim(afile)
20          format(a,' cannot open instruction file ',a,'.')
            ifail=1
            go to 9900
          endif
          read(iunit,'(a)',end=9400,err=9000) dline
          call mio_remchar(dline,achar(9))
          call utl_casetrans(dline,'lo')
          if((dline(1:3).ne.'pif').and.(dline(1:3).ne.'jif'))go to 9400
          mstruc%mrkdel(i)=dline(5:5)
          if(mstruc%mrkdel(i).eq.' ') go to 9400
50        read(iunit,'(a)',end=180,err=9000) dline
          call mio_remchar(dline,achar(9))
          if(index(dline,mstruc%mrkdel(i)).eq.0) call mio_cmprss()
          nblc=len_trim(dline)
          if(nblc.eq.0) go to 50
          if(nblc.gt.nblbmx)nblbmx=nblc
          mstruc%ninstr=mstruc%ninstr+1
          do 60 j=1,nblc
            if(dline(j:j).ne.' ') then
              if((dline(j:j).eq.'L').or.(dline(j:j).eq.'l')) mstruc%numl=mstruc%numl+1
              go to 100
            endif
60        continue
100       mstruc%asize=mstruc%asize+nblc
          go to 50
180       close(unit=iunit)
        enddo
        nblbmx=nblbmx+1
        do 300 i=1,mstruc%numoutfile
          mstruc%asize=mstruc%asize+2+len(mstruc%modoutfile(i))
300     continue
        mstruc%ninstr=mstruc%ninstr+mstruc%numoutfile
        
! -- Memory is allocated for storage of instructions.

	allocate(mstruc%a(mstruc%asize),mstruc%ll(mstruc%numl),mstruc%lcins(mstruc%ninstr),stat=ierr)
	if(ierr.ne.0)then
	  write(amessage,310)
310	  format('Error in function MIO_STORE_INSTRUCTION_SET: ',   &
          'cannot allocate sufficient memory to store instruction set.')
          ifail=1
          go to 9900
        endif
        mstruc%a=' '              ! a is an array
        
! -- The instruction set is now re-read and stored.

        ins=0
        isum=0
        do 400 i=1,mstruc%numoutfile
          call utl_addquote(mstruc%insfile(i),afile)
          iunit=utl_nextunit()
          open(unit=iunit,file=mstruc%insfile(i),status='old',iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,20) trim(errsub),trim(afile)
            ifail=1
            go to 9900
          endif
          read(iunit,*,err=9000)
          ins=ins+1
          dline(1:1)=achar(2)
          dline(2:2)=' '
          dline(3:)=mstruc%modoutfile(i)
          mstruc%lcins(ins)=isum+1
          nblc=len(mstruc%modoutfile(i))+2
          do j=1,nblc
            mstruc%a(j+isum)=dline(j:j)
          enddo
          isum=isum+nblc
350       read(iunit,322,end=181,err=9000) dline
322       format(a)
          call mio_remchar(dline,achar(9))
          if(index(dline,mstruc%mrkdel(i)).eq.0) call mio_cmprss()
          nblc=len_trim(dline)
          if(nblc.eq.0) go to 350
          ins=ins+1
          mstruc%lcins(ins)=isum+1
          do j=1,nblc
            mstruc%a(j+isum)=dline(j:j)
          enddo
          isum=isum+nblc
          go to 350
181       close(unit=iunit)
400     continue
        
        mstruc%instruction_status=1
        go to 9900

9000    write(amessage,9010) trim(errsub),trim(afile)
9010    format(a,' unable to read instruction file ',a,'.')
        ifail=1
        go to 9900
9400    write(amessage,9410) trim(errsub),trim(afile)
9410    format(a,' header of "pif" or "jif" followed by space, followed by marker delimiter ', &
        'expected on first line of instruction file ',a,'.')
        ifail=1
        go to 9900

9900    continue
        if(ifail.eq.0)then
          mio_store_instruction_set=0
        else if(ifail.eq.1)then
          mio_store_instruction_set=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return
        
end function mio_store_instruction_set



integer function mio_delete_model_output_files(estruc,mstruc)

! -- Function MIO_DELETE_MODEL_OUTPUT_FILES deletes the set of model output files pertaining
!    to the current model output interface dataset.

        implicit none

        type(err_failure_struc), intent(inout)          :: estruc
        type(mio_struc), intent(inout)                  :: mstruc

        logical                                         :: lexist
        integer                                         :: jerr,iunit,i,idir,n
        character*1                                     :: aa
        character*256                                   :: mofile

        function_name='MIO_DELETE_MODEL_OUTPUT_FILES'
        mstruc%ic_delete_model_output_files=mstruc%ic_delete_model_output_files+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

! -- Context checking is undertaken.

        if(mstruc%ic_observation_check.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_OBSERVATION_CHECK has been called.')
          ifail=1
          go to 9900
        endif

        idir=0
        if(mstruc%workdir.ne.' ') idir=1

        do i=1,mstruc%numoutfile
          if(idir.eq.0)then
            mofile=mstruc%modoutfile(i)
          else
            mofile=trim(mstruc%workdir)//trim(mstruc%modoutfile(i))
          endif
          if(utl_delete_file(mofile).ne.0)then
            call utl_addquote(mofile,afile)
            write(amessage,10) trim(afile)
10          format('Cannot delete model output file ',a,' prior to running model.')
            ifail=1
            go to 9900
          endif
        enddo

9900    continue
        if(ifail.eq.0)then
          mio_delete_model_output_files=0
        else if(ifail.eq.1)then
          mio_delete_model_output_files=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
        endif
        return

end function mio_delete_model_output_files




integer function mio_read_model_output_files(estruc,mstruc,oval)

! -- Function MIO_READ_MODEL_OUTPUT_FILES reads model output files using an instruction
!    set.

!! Important note: if an error condition occurs the error message may contain more than one line.
!! The second line will be "Instruction follows:-'. The instruction at which the error occurred
!! will follow in the string after that.

        implicit none

        type(err_failure_struc), intent(inout)             :: estruc
        type(mio_struc), intent(inout)                     :: mstruc
        double precision, intent(out), dimension(:)        :: oval

        integer                          :: ifile,il,jobs,cil,iobs,begins,ins,nblb,i,n,    &
                                            n1,n2,insnum,nblc,dumflg,marktyp,almark,iunit, &
                                            nol,jfail,insfle,mrktyp,j2,j1,n3,num1,num2,j,idir, &
                                            ilstart
        double precision                 :: rtemp
        character (len=1)                :: aa,mkrdel
        character (len=20)               :: obsnam
        character (len=15)               :: fmt
        character (len=10)               :: anum
        character (len=256)              :: flenme
        character (len=500)              :: instruction

        function_name='MIO_READ_MODEL_OUTPUT_FILES'
        mstruc%ic_read_model_output_files=mstruc%ic_read_model_output_files+1
        errsub='Error in function '//trim(function_name)//':'
        ifail=0

        instruction=' '

! -- First a context check is made.

        if(mstruc%ic_store_instruction_set.eq.0)then
          write(amessage,5) trim(errsub)
5         format(a,' this function must not be called before function MIO_STORE_INSTRUCTION_SET has been called.')
          ifail=1
          go to 9900
        endif

! -- Next the size of the observation value array is checked.

        if(size(oval).lt.mstruc%nobs)then
          write(amessage,6) trim(errsub)
6         format(a,' size of supplied OVAL array is too small.')
          ifail=1
          go to 9900
        endif

! -- Now processing commences.

        oval=-1.1d270    ! an array   MD...Scientific notation for double precision assignment

        errsub='Error reading model output file(s):'
        ifail=0
        idir=0
        if(mstruc%workdir.ne.' ') idir=1

        mstruc%mcall=mstruc%mcall+1
        ifile=0
        il=0
        jobs=0
        mkrdel=mstruc%mrkdel(1)
        cil=0
        iobs=1
        begins=0

        ins=1
10      if(ins.lt.mstruc%ninstr)then
          nblb=mstruc%lcins(ins+1)-mstruc%lcins(ins)
        else
          nblb=mstruc%asize-mstruc%lcins(ins)+1
        endif
        instruction=' '
        do 20 i=1,nblb
          instruction(i:i)=mstruc%a(mstruc%lcins(ins)+i-1)
20      continue
25      n2=0
        insnum=0

50      call mio_getint(jfail,instruction,n1,n2,nblb,mkrdel)
        if(jfail.ne.0)then
          write(amessage,49) trim(errsub)
49        format(a,' missing marker delimiter in user-supplied instruction.')
          go to 9995
        endif
51      if(n1.eq.0) go to 1000
        insnum=insnum+1
        if(insnum.eq.1)then
          if(instruction(n1:n1).ne.'&') then
            mrktyp=0
            almark=1
            begins=0
            ilstart=il
          else
            if(ins.eq.insfle+1) then
              write(amessage,52) trim(errsub)
52            format(a,' first instruction line in instruction file cannot start ', &
              'with continuation character.')
              go to 9995
            endif
            if(begins.eq.1)then
              ins=ins-1
              go to 10
            endif
          endif
        endif
        if(ichar(instruction(n1:n1)).eq.2)then
          if(ifile.ne.0) close(unit=iunit)
          do 60 i=n1+1,nblb
            if(instruction(i:i).ne.' ') go to 70
60        continue
70        flenme=instruction(i:nblb)
          if(idir.eq.1) flenme=trim(mstruc%workdir)//trim(flenme)
          iunit=utl_nextunit()
          do i=1,4
            open(unit=iunit,file=flenme,status='old',iostat=ierr)
            if(ierr.eq.0) exit
            call utl_wait_std(100)
          enddo
          call utl_addquote(flenme,afile)
          if(ierr.ne.0)then
            write(amessage,71) trim (errsub),trim(afile)
71          format(a,' cannot open model output file ',a,'.')
            instruction=' '
            go to 9995
          endif
          ifile=ifile+1
          cil=0
          mkrdel=mstruc%mrkdel(ifile)
          insfle=ins
          go to 1000
        else if((instruction(n1:n1).eq.'l').or.(instruction(n1:n1).eq.'L'))then
          if(il.ne.ilstart)then
            write(amessage,72) trim(errsub)
72          format(a,' line advance item can only occur at the beginning of an instruction line.')
            go to 9995
          endif
          almark=0
          il=il+1
          if(mstruc%mcall.eq.1)then
            if(n2.le.n1) go to 9050     ! put in pest
            write(fmt,150) n2-n1
150         format('(i',i4,')')
            read(instruction(n1+1:n2),fmt,err=9050) nol
            mstruc%ll(il)=nol
          else
            nol=mstruc%ll(il)
          endif
          if(nol.gt.1) then
            do 160 i=1,nol-1
              read(iunit,*,end=9100,err=9850)
              cil=cil+1
160         continue
          endif
          read(iunit,22,end=9100,err=9850) dline
22        format(a)
          if(index(dline,char(9)).ne.0) call mio_tabrep()
          cil=cil+1
          nblc=len_trim(dline)
          mrktyp=1
          j1=0
        else if(instruction(n1:n1).eq.mkrdel)then
          if(mrktyp.eq.0)then
200         read(iunit,22,end=9100,err=9850) dline
            if(index(dline,char(9)).ne.0) call mio_tabrep()
            cil=cil+1
            j1=index(dline,instruction(n1+1:n2-1))
            if(j1.eq.0) go to 200
            nblc=len_trim(dline)
            j1=j1+n2-n1-2
            mrktyp=1
          else
            if(j1.ge.nblc) then
              if(almark.eq.1) then
                begins=1
                go to 25
              endif
              go to 9200
            endif
            j2=index(dline(j1+1:nblc),instruction(n1+1:n2-1))
            if(j2.eq.0) then
              if(almark.eq.1) then
                begins=1
                go to 25
              endif
              go to 9200
            endif
            j1=j1+j2
            j1=j1+n2-n1-2
          endif
        else if(instruction(n1:n1).eq.'&')then
          if(insnum.ne.1) then
            write(amessage,201) trim(errsub)
201         format(a,' if present, continuation character must be first instruction on ', &
            'an instruction line.')
            go to 9995
          endif
        else if((instruction(n1:n1).EQ.'w').or.(instruction(n1:n1).eq.'W'))then
          almark=0
          if(j1.ge.nblc) go to 9400
          j2=index(dline(j1+1:nblc),' ')
          if(j2.eq.0) go to 9400
          j1=j1+j2
          do 210 i=j1,nblc
            if(dline(i:i).ne.' ') go to 220
210       continue
          i=nblc+1
220       j1=i-1
        else if((instruction(n1:n1).eq.'t').or.(instruction(n1:n1).eq.'T'))then
          almark=0
          if(n2.le.n1) go to 9000       ! put in PEST
          write(fmt,150) n2-n1
          read(instruction(n1+1:n2),fmt,err=9000) j2
          if(j2.lt.j1) then
            call utl_int2char(cil,anum)
            write(amessage,221) trim(errsub),trim(anum),trim(afile)
221         format(a,' backwards move to tab position not allowed on line ',a,  &
            ' of model output file ',a,'.')
            go to 9995
          endif
          j1=j2
          if(j1.gt.nblc) then
            call utl_int2char(cil,anum)
            write(amessage,222) trim(errsub),trim(anum),trim(afile)
222         format(a,' tab position beyond end of line at line ',a,' of ', &
            'model output file ',a,'.')
            go to 9995
          endif
        else if((instruction(n1:n1).eq.'[').or.(instruction(n1:n1).eq.'('))then
          almark=0
          aa=instruction(n1:n1)
          jobs=jobs+1
          if(mstruc%mcall.eq.1)then
            if(aa.eq.'[')then
              n3=index(instruction(n1:n2),']')
            else
              n3=index(instruction(n1:n2),')')
            endif
            if(n3.eq.0)then
              call utl_int2char(cil,anum)
              write(amessage,226) trim(errsub)
226           format(a,' missing "]" or ")" character in instruction.')
              go to 9995
            endif
            n3=n3+n1-1
            obsnam=instruction(n1+1:n3-1)
            call mio_which1(jfail,mstruc%nobs,iobs,mstruc%aobs,obsnam)
            if(jfail.ne.0) go to 9700
            call mio_getnum(jfail,instruction,n3,n2,num1,num2,fmt)
            IF(jfail.ne.0) then
              write(amessage,223) trim(errsub)
223           format(a,' cannot interpret user-supplied instruction for reading model ', &
              'output file.')
              go to 9995
            endif
            mstruc%obsn1(jobs)=num1
            mstruc%obsn2(jobs)=num2
            mstruc%iiobs(jobs)=iobs
          else
            num1=mstruc%obsn1(jobs)
            num2=mstruc%obsn2(jobs)
            iobs=mstruc%iiobs(jobs)
          endif
          if(aa.eq.'(') then
            call mio_gettot(jfail,num1,num2,nblc)
            if(jfail.ne.0)then
              call utl_int2char(cil,anum)
              write(amessage,224) trim (errsub),trim(mstruc%aobs(iobs)),trim(anum),   &
              trim(afile)
224           format(a,' cannot find observation "',a,'" on line ',a,     &
              ' of model output file ',a,'.')
              go to 9995
            endif
          else
            if(num1.gt.nblc)then
              call utl_int2char(cil,anum)
              write(amessage,224) trim(errsub),trim(mstruc%aobs(iobs)),trim(anum),trim(afile)
              go to 9995
            endif
            if(num2.gt.nblc) num2=nblc
            if(dline(num1:num2).eq.' ')then
              call utl_int2char(cil,anum)
              write(amessage,224) trim(errsub),trim(mstruc%aobs(iobs)),trim(anum),trim(afile)
              go to 9995
            endif
          endif
          write(fmt,250) num2-num1+1
250       format('(f',i4,'.0)')
          if(oval(iobs).gt.-1.0d270) go to 9870
          read(dline(num1:num2),fmt,err=260) oval(iobs)
          j1=num2
          go to 50
260       continue
          call utl_int2char(cil,anum)
          write(amessage,261) trim(errsub),trim(mstruc%aobs(iobs)),trim(anum),trim(afile)
261       format(a,' cannot read observation "',a,'" from line ',a,     &
          ' of model output file ',a,'.')
          go to 9995
        else if(instruction(n1:n1).eq.'!') then
          almark=0
          call utl_casetrans(instruction(n1+1:n2-1),'lo')
          if((n2-n1.ne.4).or.(instruction(n1+1:n2-1).ne.'dum'))then
            jobs=jobs+1
            if(mstruc%mcall.eq.1) then
              obsnam=instruction(n1+1:n2-1)
              call mio_which1(jfail,mstruc%nobs,iobs,mstruc%aobs,obsnam)
              if(jfail.ne.0) go to 9700
              mstruc%iiobs(jobs)=iobs
            else
              iobs=mstruc%iiobs(jobs)
            endif
            dumflg=0
          else
            dumflg=1
          endif
          call mio_getnxt(jfail,j1,num1,num2,nblc)
          if(jfail.ne.0) then
            if(dumflg.eq.0) then
              call utl_int2char(cil,anum)
              write(amessage,224) trim(errsub),trim(mstruc%aobs(iobs)),trim(anum),trim(afile)
              go to 9995
            else
              call utl_int2char(cil,anum)
              write(amessage,224) trim(errsub),'dum',trim(anum),trim(afile)
              go to 9995
            endif
          endif
          write(fmt,250) num2-num1+1
          read(dline(num1:num2),fmt,err=270) rtemp
          if(dumflg.eq.0)then
            if(oval(iobs).gt.-1.0d270) go to 9870
            oval(iobs)=rtemp
          endif
          j1=num2
          go to 50
270       call mio_getint(jfail,instruction,n1,n2,nblb,mkrdel)
          if(jfail.ne.0) then
            write(amessage,271) trim(errsub)
271         format(a,' missing marker delimiter in user-supplied instruction set.')
            go to 9995
          endif  
          if(n1.eq.0)then
            if(dumflg.eq.1) go to 9950
            go to 9800
          endif
          if(instruction(n1:n1).ne.mkrdel) then
            if(dumflg.eq.1) go to 9950
            go to 9800
          endif
          j2=index(dline(j1+1:nblc),instruction(n1+1:n2-1))
          if(j2.eq.0) then
            if(dumflg.eq.1) go to 9950
            go to 9800
          endif
          num2=j1+j2-1
          if(num2.lt.num1)then
            if(dumflg.eq.1) go to 9950
            go to 9800
          endif
          write(fmt,250) num2-num1+1
          if(dumflg.eq.1)then
            read(dline(num1:num2),fmt,err=9950) rtemp
          else
            if(oval(iobs).gt.-1.0d270) go to 9870
            read(dline(num1:num2),fmt,err=9800) oval(iobs)
          endif
          j1=num2
          go to 51
        else
          write(amessage,272) trim(errsub)
272       format(a,' cannot interpret user-supplied instruction for reading model ',  &
          'output file.')
          go to 9995
        endif
        go to 50
1000    ins=ins+1
        if(ins.le.mstruc%ninstr) go to 10

        if(mstruc%mcall.eq.1)then
          do 1100 i=1,mstruc%nobs
          do 1050 j=1,jobs
          if(mstruc%iiobs(j).eq.i) go to 1100
1050      continue
          write(amessage,1051) trim(errsub),trim(mstruc%aobs(i))
1051      format(a,' observation "',a,'" not referenced in the user-supplied instruction set.')
          instruction=' '
          go to 9995
1100      continue
        endif

        close(unit=iunit)
        instruction=' '

        go to 9900

9000    write(amessage,9010) trim(errsub)
9010    format(a,' cannot read tab position from user-supplied instruction.')
        go to 9995
9050    write(amessage,9060) trim(errsub)
9060    format(a,' cannot read line advance item from user-supplied instruction.')
        go to 9995
9100    write(amessage,9110) trim(errsub),trim(afile)
9110    format(a,' unexpected end to model output file ',a,'.')
        go to 9995
9200    call utl_int2char(cil,anum)
        write(amessage,9210) trim(errsub),trim(anum),trim(afile)
9210    format(a,' unable to find secondary marker on line ',a,   &
        ' of model output file ',a,'.')
        go to 9995
9400    call utl_int2char(cil,anum)
        write(amessage,9410) trim(errsub),trim(anum),trim(afile)
9410    format(a,' unable to find requested whitespace, or whitespace ',  &
        'precedes end of line at line ',a,' of model output file ',a,'.')
        go to 9995
9700    write(amessage,9710) trim(errsub),trim(obsnam)
9710    format(a,' observation name "',a,'" from user-supplied instruction set ',&
        'is not cited in instruction list.')
        go to 9995
9800    call utl_int2char(cil,anum)
        write(amessage,9810) trim(errsub),trim(mstruc%aobs(iobs)),trim(anum),trim(afile)
9810    format(a,' cannot read observation "',a,'" from line ',a,   &
        ' of model output file ',a,'.')
!        instruction=' '
        go to 9995
9850    write(amessage,9860) trim(afile)
9860    format('Unable to read model output file ',a,'.')
        instruction=' '
        go to 9995
9870    write(amessage,9880) trim(errsub),trim(mstruc%aobs(iobs))
9880    format(a,' observation "',a,'" already cited in instruction set.')
        go to 9995
9950    call utl_int2char(cil,anum)
        write(amessage,9810) trim(errsub),'dum',trim(anum),trim(afile)
!        instruction=' '
        go to 9995

9995    ifail=1
        mstruc%mcall=mstruc%mcall-1
        if(mstruc%mcall.lt.0) mstruc%mcall=0
        go to 9900

9900    continue
        if(ifail.eq.0)then
          mio_read_model_output_files=0
        else if(ifail.eq.1)then
          mio_read_model_output_files=1
          call err_reset(estruc)
          call err_add_error(estruc,amessage,function_name)
          if(instruction.ne.' ')then
            call err_add_error(estruc,'Instruction follows:-',function_name)
            call err_add_error(estruc,instruction,function_name)
          endif
        endif
        return

end function mio_read_model_output_files



integer function mio_finalise(estruc,mstruc)

! -- Function MIO_FINALISE de-allocates memory usage by the model_input_output
!    module.

        implicit none
      
        type(err_failure_struc), intent(inout)             :: estruc
        type(mio_struc), intent(inout)                     :: mstruc

        function_name='MIO_FINALISE'
        ifail=0
        mstruc%ic_finalise=mstruc%ic_finalise+1
        mstruc%ic_initialise=0

        mstruc%template_status=0
        mstruc%instruction_status=0
      
        if(associated(mstruc%tempfile))deallocate(mstruc%tempfile,stat=ierr)
        if(ierr.ne.0) then
          ifail=1
        else
          nullify(mstruc%tempfile)
        endif
        if(associated(mstruc%modinfile))deallocate(mstruc%modinfile,stat=ierr)
        if(ierr.ne.0) then
          ifail=1
        else
          nullify(mstruc%modinfile)
        endif
        if(associated(mstruc%pardel))deallocate(mstruc%pardel,stat=ierr)
        if(ierr.ne.0) then
          ifail=1
        else
          nullify(mstruc%pardel)
        endif
        if(associated(mstruc%insfile))deallocate(mstruc%insfile,stat=ierr)
        if(ierr.ne.0) then
          ifail=1
        else
          nullify(mstruc%insfile)
        endif
        if(associated(mstruc%modoutfile))deallocate(mstruc%modoutfile,stat=ierr)
        if(ierr.ne.0) then
          ifail=1
        else
          nullify(mstruc%modoutfile)
        endif
        if(associated(mstruc%mrkdel))deallocate(mstruc%mrkdel,stat=ierr)
        if(ierr.ne.0) then
          ifail=1
        else
          nullify(mstruc%mrkdel)
        endif
        if(associated(mstruc%nw))deallocate(mstruc%nw,stat=ierr)
        if(ierr.ne.0) then
          ifail=1
        else
          nullify(mstruc%nw)
        endif
        if(associated(mstruc%pword))deallocate(mstruc%pword,stat=ierr)
        if(ierr.ne.0) then
          ifail=1
        else
          nullify(mstruc%pword)
        endif
        if(associated(mstruc%obsn1))deallocate(mstruc%obsn1,stat=ierr)
        if(ierr.ne.0) then
          ifail=1
        else
          nullify(mstruc%obsn1)
        endif
        if(associated(mstruc%obsn2))deallocate(mstruc%obsn2,stat=ierr)
        if(ierr.ne.0) then
          ifail=1
        else
          nullify(mstruc%obsn2)
        endif
        if(associated(mstruc%iiobs))deallocate(mstruc%iiobs,stat=ierr)
        if(ierr.ne.0) then
          ifail=1
        else
          nullify(mstruc%iiobs)
        endif
        if(associated(mstruc%a))deallocate(mstruc%a,stat=ierr)
        if(ierr.ne.0) then
          ifail=1
        else
          nullify(mstruc%a)
        endif
        if(associated(mstruc%ll))deallocate(mstruc%ll,stat=ierr)
        if(ierr.ne.0) then
          ifail=1
        else
          nullify(mstruc%ll)
        endif
        if(associated(mstruc%lcins))deallocate(mstruc%lcins,stat=ierr)
        if(ierr.ne.0) then
          ifail=1
        else
          nullify(mstruc%lcins)
        endif
        if(associated(mstruc%apar))deallocate(mstruc%apar,stat=ierr)
        if(ierr.ne.0) then
          ifail=1
        else
          nullify(mstruc%apar)
        endif
        if(associated(mstruc%aobs))deallocate(mstruc%aobs,stat=ierr)
        if(ierr.ne.0) then
          ifail=1
        else
          nullify(mstruc%aobs)
        endif

9900    continue
        if(ifail.eq.0)then
          mio_finalise=0
        else if(ifail.eq.1)then
          mio_finalise=0
        endif

        return

end function mio_finalise







subroutine mio_getnxt(ifail,j1,num1,num2,nblc)

! -- Subroutine MIO_GETNXT gets the next space-delimited word.

        implicit none        

        integer ifail
        integer j1,num1,num2,nblc,i

        ifail=0
        do 20 i=j1+1,nblc
        if(dline(i:i).ne.' ') go to 50
20      continue
        ifail=1
        return
50      num1=i
        i=index(dline(num1:nblc),' ')
        if(i.eq.0) then
          num2=nblc
        else
          num2=num1+i-2
        endif

        return

end subroutine mio_getnxt


subroutine mio_gettot(ifail,j1,j2,nblc)

! -- Subroutine MIO_GETTOT determines the exact position occupied by a number.

        implicit none
        integer ifail
        integer j1,j2,nblc,i
 
        ifail=0
        if(j1.gt.nblc)then
          ifail=1
          return
        endif
        if(j2.gt.nblc)j2=nblc
        if(dline(j2:j2).eq.' ') then
          do 10 i=j2,j1,-1
          if(dline(i:i).ne.' ')then
            j2=i
            go to 100
          endif
10        continue
          ifail=1
          return
        else
          if(j2.eq.nblc) go to 100
          do 20 i=j2,nblc
          if(dline(i:i).eq.' ') then
            j2=i-1
            go to 100
          endif
20        continue
          j2=nblc
        endif
100     if(j1.eq.1) go to 200
        do 120 i=j1,1,-1
        if(dline(i:i).eq.' ') then
          j1=i+1
          go to 200
        endif
120     continue
        j1=1
200     return

end subroutine mio_gettot


subroutine mio_getint(ifail,buf,n1,n2,nblb,mrkdel)

! -- Subroutine MIO_GETINT gets the next stored instruction for processing.

        integer n1,n2,nblb,i,ii
        integer ifail
        character mrkdel
        character*(*) buf

        ifail=0
        if(n2.ge.nblb) then
          n1=0
          return
        endif
        do 10 i=n2+1,nblb
        if((buf(i:i).ne.' ').and.(ichar(buf(i:i)).ne.9)) go to 50
10      continue
        n1=0
        return
50      n1=i
        if(buf(n1:n1).ne.mrkdel)then
          i=index(buf(n1:nblb),' ')
          ii=index(buf(n1:nblb),char(9))
          if((i.eq.0).and.(ii.eq.0))then
            i=0
          else if(i.eq.0)then
            i=ii
          else if(ii.eq.0) then
            i=i
          else
            i=min(i,ii)
          endif
          if(i.ne.0) then
            n2=n1+i-2
          else
            n2=nblb
          endif
        else
          if(n1.eq.nblb)then
            ifail=1
            return
          endif
          i=index(buf(n1+1:nblb),mrkdel)
          if(i.eq.0) then
            ifail=1
            return
          endif
          n2=n1+i
        endif

        return
 end subroutine mio_getint


subroutine mio_tabrep()

! -- Subroutine MIO_TABREP replaces a tab by blank space(s).

        integer llen,i,j,k,nblc

        llen=len(dline)
        do 10 i=llen,1,-1
        if(dline(i:i).ne.' ') go to 20
10      continue
        return
20      nblc=i

        i=0
30      i=i+1
        if(i.gt.nblc)return
        if(ichar(dline(i:i)).ne.9) go to 30
        j=((i-1)/8+1)*8-i
        if(j.eq.0) then
          dline(i:i)=' '
        else
          dline(i:i)=' '
          nblc=nblc+j
          if(nblc.gt.llen) nblc=llen
          do 50 k=nblc,((i-1)/8+1)*8,-1
          dline(k:k)=dline(k-j:k-j)
50        continue
          do 60 k=i+1,min(nblc,i+j)
          dline(k:k)=' '
60        continue
          i=i+j
        endif
        go to 30

end subroutine mio_tabrep



subroutine mio_parnam(ifail,j1,j2,tpar)

! -- Subroutine MIO_PARNAM extracts a parameter name from a string.

        implicit none
        integer, intent(out)          :: ifail   ! reports error condition
        integer, intent(in)           :: j1,j2   ! beginning and end of word
        character (len=*), intent(out):: tpar    ! the extracted parameter
        
        integer             :: i,j

        ifail=0
        tpar=' '
        if(j2-j1.le.1) then
          ifail=1
          return
        endif
        do 10 i=j1+1,j2-1
        if(dline(i:i).eq.' ') go to 10
        go to 30
10      continue
        ifail=2
        return
30      j=min(12,j2-i)
        tpar(1:j)=dline(i:i+j-1)
        return

end subroutine mio_parnam






subroutine mio_cmprss()

! -- Subroutine MIO_CMPRSS compresses an instruction line by removing excess
! -- blank characters.

        implicit none
        
        integer nblc,j

        if(dline.eq.' ') return
10      nblc=len_trim(dline)
        j=index(dline(1:nblc),'  ')
        if(j.ne.0) then
          dline(j+1:)=adjustl(dline(j+1:))
          go to 10
        endif
        return

end subroutine mio_cmprss


subroutine mio_getnum(ifail,buf,n3,n2,num1,num2,fmt)

! -- Subroutine MIO_GETNUM retrieves character positions from fixed and
!    semi-fixed observation instructions.

        integer n3,num1,num2,i,n2
        integer ifail
        character*(*) buf
        character*(*) fmt

        ifail=0
        i=index(buf(n3+1:n2),':')
        if(i.eq.0) go to 100
        write(fmt,20) i-1
20      format('(i',i3,')')
        read(buf(n3+1:n3+i-1),fmt,err=100) num1
        n3=n3+i
        i=n2-n3
        if(i.lt.1) go to 100
        write(fmt,20) i
        read(buf(n3+1:n2),fmt,err=100) num2
        return
100     ifail=1
        return

end subroutine mio_getnum




subroutine mio_which1(ifail,ndim,idim,aname,tname)

! -- Subroutine MIO_WHICH1 finds a string in an array of strings.

        implicit none

        integer, intent(out)                            :: ifail    ! error indicator
        integer, intent(in)                             :: ndim     ! number of elements
        integer, intent(inout)                          :: idim     ! where to start the search
        character (len=*), intent(in), dimension(ndim)  :: aname    ! an array of object names
        character (len=*), intent(inout)                :: tname    ! object name to look for

        integer                             :: i

        ifail=0
        if((idim.lt.1).or.(idim.gt.ndim)) idim=1
        call utl_casetrans(tname,'lo')
        if(tname.eq.aname(idim)) return
        if(idim.ne.ndim)then
          do 20 i=idim+1,ndim
            if(tname.eq.aname(i))then
              idim=i
              return
            endif
20        continue
        endif
        if(idim.ne.1)then
          do 40 i=idim-1,1,-1
          if(tname.eq.aname(i)) then
            idim=i
            return
          endif
40        continue
        endif
        ifail=1
        return

end subroutine mio_which1




subroutine mio_remchar(astring,ach)

       implicit none

       character*(*), intent(inout) :: astring
       character*(*), intent(in)    :: ach

       integer ll,ii,icount

       icount=0
       ll=len_trim(ach)

10     ii=index(astring,ach)
       if(ii.eq.0) then
         if(icount.eq.0)return
         go to 20
       endif
       icount=icount+1
       astring(ii:ii-1+ll)=' '
       go to 10

20     astring=adjustl(astring)
       return

end subroutine mio_remchar



subroutine mio_wrtsig(ifail,val,word,nw,precis,tval,nopnt)
! --
! -- Subroutine WRTSIG writes a number into a confined space with maximum
! -- precision.
! --

!       failure criteria:
!           ifail= 1 ...... number too large or small for single precision type
!           ifail= 2 ...... number too large or small for double precision type
!           ifail= 3 ...... field width too small to represent number
!           ifail=-1 ...... internal error type 1
!           ifail=-2 ...... internal error type 2
!           ifail=-3 ...... internal error type 3

        integer precis,lw,pos,inc,d,p,w,j,jj,k,jexp,n,jfail,nw,epos,pp,nopnt,kexp,iflag,lexp
        integer ifail
        double precision val,tval
        character*29 tword,ttword,fmt*14
        character*(*) word

!       The following line overcomes what appears to be a bug in the LF90
!       compiler

#ifdef LAHEY
        if(abs(val).lt.1.0d-300) val=0.0d0
#endif

        lexp=0
        iflag=0
        word=' '
        pos=1
        if(val.lt.0.0d0)pos=0
#ifdef USE_D_FORMAT
        write(tword,'(1PD23.15D3)') val
#else
        write(tword,'(1PE23.15E3)') val
#endif
        call utl_casetrans(tword,'hi')
        read(tword(20:23),'(i4)') jexp
        epos=1
        if(jexp.lt.0)epos=0

        jfail=0
        ifail=0
        if(precis.eq.0)then
          lw=min(15,nw)
        else
          lw=min(23,nw)
        endif

        n=0
        if(nopnt.eq.1)n=n+1
        if(pos.eq.1)n=n+1
        if(precis.eq.0)then
          if(abs(jexp).gt.38)then
            ifail=1
            return
          endif
          if(pos.eq.1) then
            if(lw.ge.13) then
              write(word,'(1pe14.7)',err=80) val
              go to 200
            endif
          else
            if(lw.ge.14)then
              write(word,'(1pe14.7)',err=80) val
              go to 200
            endif
          endif
          if(lw.ge.14-n) then
            lw=14-n
            go to 80
          endif
        else
          if(abs(jexp).gt.275)then
            ifail=2
            return
          endif
          if(pos.eq.1) then
            if(lw.ge.22) then
#ifdef USE_D_FORMAT
              write(word,'(1PD22.15D3)',err=80) val
#else
              write(word,'(1PE22.15E3)',err=80) val
#endif
              go to 200
            endif
          else
            if(lw.ge.23) then
#ifdef USE_D_FORMAT
              write(word,'(1PD23.15D3)',err=80) val
#else
              write(word,'(1PE23.15E3)',err=80) val
#endif
              go to 200
            endif
          endif
          if(lw.ge.23-n)then
            lw=23-n
            go to 80
          endif
        endif

        if(nopnt.eq.1)then
          if((jexp.eq.lw-2+pos).or.(jexp.eq.lw-3+pos))then
            write(fmt,15)lw+1
15          format('(f',i2,'.0)')
            write(word,fmt,err=19) val
            if(index(word,'*').ne.0) go to 19
            if(word(1:1).eq.' ') go to 19
            word(lw+1:lw+1)=' '
            go to 200
          endif
        endif
19      d=min(lw-2+pos,lw-jexp-3+pos)
20      if(d.lt.0) go to 80
        write(fmt,30) lw,d
30      format('(f',i2,'.',i2,')')
        write(word,fmt,err=80) val
        if(index(word,'*').ne.0) then
          d=d-1
          go to 20
        endif
        k=index(word,'.')
        if(k.eq.0)then
          ifail=-1
          return
        endif
        if((k.eq.1).or.((pos.eq.0).and.(k.eq.2)))then
          do 70 j=1,3
          if(k+j.gt.lw) go to 75
          if(word(k+j:k+j).ne.'0') go to 200
70        continue
          go to 80
75        ifail=3
          return
        endif
        go to 200

80      word=' '
        if(nopnt.eq.0)then
          d=lw-7
          if(pos.eq.1) d=d+1
          if(epos.eq.1) d=d+1
          if(abs(jexp).lt.100) d=d+1
          if(abs(jexp).lt.10) d=d+1
          if((jexp.ge.100).and.(jexp-(d-1).lt.100))then
            p=1+(jexp-99)
            d=d+1
            lexp=99
          else if((jexp.ge.10).and.(jexp-(d-1).lt.10))then
            p=1+(jexp-9)
            d=d+1
            lexp=9
          else if((jexp.eq.-10).or.(jexp.eq.-100)) then
            iflag=1
            d=d+1
          else
            p=1
          endif
          inc=0
85        if(d.le.0) go to 300
          if(iflag.eq.0)then
            write(fmt,100,err=300) p,d+7,d-1
          else
            write(fmt,100,err=300) 0,d+8,d
          endif
          write(tword,fmt) val
          call utl_casetrans(tword,'hi')
          if(iflag.eq.1) go to 87
          read(tword(d+4:d+7),'(i4)',err=500) kexp
          if(((kexp.eq.10).and.((jexp.eq.9).or.(lexp.eq.9))).or.     &
          ((kexp.eq.100).and.((jexp.eq.99).or.lexp.eq.99))) then
            if(inc.eq.0)then
              if(lexp.eq.0)then
                if(d-1.eq.0) then
                  d=d-1
                else
                  p=p+1
                endif
              else if(lexp.eq.9)then
                if(jexp-(d-2).lt.10) then
                  p=p+1
                else
                  d=d-1
                endif
              else if(lexp.eq.99)then
                if(jexp-(d-2).lt.100)then
                  p=p+1
                else
                  d=d-1
                endif
              endif
              inc=inc+1
              go to 85
            endif
          endif
#ifdef USE_D_FORMAT
87        j=index(tword,'D')
#else
87        j=index(tword,'E')
#endif
          go to 151
        endif
        inc=0
        p=lw-2
        pp=jexp-(p-1)
        if(pp.ge.10)then
          p=p-1
          if(pp.ge.100)p=p-1
        else if(pp.lt.0)then
          p=p-1
          if(pp.le.-10)then
            p=p-1
            if(pp.le.-100)p=p-1
          endif
        endif
        if(pos.eq.0)p=p-1
90      continue
        d=p-1
        w=d+8
        write(fmt,100) p,w,d
        if(d.lt.0)then
          if(jfail.eq.1) go to 300
          jfail=1
          p=p+1
          go to 90
        endif
#ifdef USE_D_FORMAT
100     format('(',I2,'pD',I2,'.',I2,'D3)')
#else
100     format('(',I2,'pE',I2,'.',I2,'E3)')
#endif
        write(tword,fmt) val
        call utl_casetrans(tword,'hi')
#ifdef USE_D_FORMAT
        j=index(tword,'D')
#else
        j=index(tword,'E')
#endif
        if(tword(j-1:j-1).ne.'.')then
          ifail=-1
          return
        endif
        n=1
        if(tword(j+1:j+1).eq.'-') n=n+1
        if(tword(j+2:j+2).ne.'0') then
          n=n+2
          go to 120
        endif
        if(tword(j+3:j+3).ne.'0') n=n+1
120     n=n+1
        if(j+n-2-pos.lt.lw)then
          if(inc.eq.-1) go to 150
          ttword=tword
          p=p+1
          inc=1
          go to 90
        else if(j+n-2-pos.eq.lw) then
          go to 150
        else
          if(inc.eq.1)then
            tword=ttword
            go to 150
          endif
          if(jfail.eq.1) go to 300
          p=p-1
          inc=-1
          go to 90
        endif

150     j=index(tword,'.')
151     if(pos.eq.0)then
          k=1
        else
         k=2
        endif
        word(1:j-k)=tword(k:j-1)
        jj=j
        j=j-k+1
        if(precis.eq.0)then
          word(j:j)='E'
        else
          word(j:j)='D'
        endif
        jj=jj+2
        if(nopnt.eq.0) jj=jj-1
        if(tword(jj:jj).eq.'-')then
          j=j+1
          word(j:j)='-'
        endif
        if(tword(jj+1:jj+1).ne.'0')then
          j=j+2
          word(j-1:j)=tword(jj+1:jj+2)
          go to 180
        endif
        if(tword(jj+2:jj+2).ne.'0')then
          j=j+1
          word(j:j)=tword(jj+2:jj+2)
        endif
180     j=j+1
        word(j:j)=tword(jj+3:jj+3)
        if(iflag.eq.1)then
          if(pos.eq.1)then
            jj=1
          else
            jj=2
          endif
          n=len_trim(word)
          do 190 j=jj,n-1
190       word(j:j)=word(j+1:j+1)
          word(n:n)=' '
        endif

200     if(len_trim(word).gt.lw)then
          ifail=-2
          return
        endif
        write(fmt,30) lw,0
        read(word,fmt,err=400) tval
        return
300     ifail=3
        return
400     ifail=-3
        return
500     ifail=-2
        return

end subroutine mio_wrtsig


end module model_input_output





! Notes:-

! -- As presently programmed, this module allows a parameter to be unrepresented in all model input
!      template files. However all observations must be cited in instruction files. This can be easily
!      altered of course.

