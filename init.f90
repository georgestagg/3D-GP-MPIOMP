module init
	use workspace
	use io
	contains
	subroutine initCond
		implicit none
		integer :: f
		if(RANK .eq. 0) then
			write (6, *) 'Applying initial condition...'
		end if
		TIME = 0.0d0
		if(initialCondType .eq. 0) then
			do f = 1,FLUIDS
				WS%FLUID(f)%GRID = 1.0d0
			end do
		else if (initialCondType .eq. 1) then
			do f = 1,FLUIDS
				call makeRandomPhase(WS%FLUID(f))
			end do
		else if (initialCondType .eq. 2) then
			call makeNonEquibPhase
		else if (initialCondType .eq. 3) then
			call read_wf_file(ICRfilename)
		else
			do f = 1,FLUIDS
				WS%FLUID(f)%GRID = 1.0d0
			end do
		end if
	end subroutine
end module