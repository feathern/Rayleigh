!
!  Copyright (C) 2018 by the authors of the RAYLEIGH code.
!
!  This file is part of RAYLEIGH.
!
!  RAYLEIGH is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 3, or (at your option)
!  any later version.
!
!  RAYLEIGH is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with RAYLEIGH; see the file LICENSE.  If not see
!  <http://www.gnu.org/licenses/>.
!

Module Messages
    Use Buffered_Output, Only : stdout

    Implicit None
Contains
    Subroutine Initialize_Messages
            Call init_random_seed()
    End Subroutine Initialize_Messages

    Subroutine Random_Integer(m,n)
            Real :: r
            Integer, Intent(Out) :: n
            Integer, Intent(In) :: m
            Call RANDOM_NUMBER(r)
            n = FLOOR((m+1)*r) 
    End Subroutine Random_Integer

    Subroutine init_random_seed()
        ! Provided by gnu.org
        Integer :: i, n, clock
        Integer, Dimension(:), Allocatable :: seed

        Call RANDOM_SEED(size = n)
        Allocate(seed(n))

        Call SYSTEM_CLOCK(count=clock)

        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        Call RANDOM_SEED(put=seed)

        Deallocate(seed)
    END Subroutine

    Subroutine Encouraging_Message(msg)
        Intent(Out) :: msg
        Integer :: msg_n
        msg_n = 
    End Subroutine Random_Encouraging_Message
End Module Messages
