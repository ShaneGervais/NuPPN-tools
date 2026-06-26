program update_netgen_rate_from_dat
    use iso_fortran_env, only: real64, error_unit
    implicit none

    integer, parameter :: line_len = 512
    integer, parameter :: token_len = 128
    integer, parameter :: label_len = 96
    integer, parameter :: max_tokens = 128
    integer, parameter :: max_points = 4096

    character(len=line_len), allocatable :: lines(:)
    character(len=line_len) :: dat_path, netgen_path, reaction_arg, out_path, requested_out_path
    character(len=label_len) :: reaction_key
    real(real64), allocatable :: source_t9(:), source_rates(:)
    integer :: nargs, nlines, nsource, rate_start, rate_end

    nargs = command_argument_count()
    if (nargs < 3 .or. nargs > 4) then
        call print_usage()
        stop 2
    end if

    call get_command_argument(1, dat_path)
    call get_command_argument(2, netgen_path)
    call get_command_argument(3, reaction_arg)
    if (nargs == 4) then
        call get_command_argument(4, requested_out_path)
    else
        requested_out_path = ''
    end if

    reaction_key = canonical_key(trim(reaction_arg))
    out_path = normalized_output_path(trim(requested_out_path), trim(netgen_path))
    call ensure_project_tmp_dir(out_path)

    call read_article_rates(trim(dat_path), reaction_key, source_t9, source_rates, nsource)
    call read_text_file(trim(netgen_path), lines, nlines)
    call find_netgen_rate_block(lines, nlines, reaction_key, rate_start, rate_end)
    call replace_rate_block(lines, rate_start, rate_end, source_t9, source_rates, nsource)
    call write_text_file(trim(out_path), lines, nlines)

    write(*,'(A)') 'updated reaction: ' // trim(reaction_arg)
    write(*,'(A,I0,A,I0)') 'replaced netgen lines: ', rate_start, '-', rate_end
    write(*,'(A,I0)') 'article rate points used: ', nsource
    write(*,'(A)') 'wrote: ' // trim(out_path)

contains

    subroutine print_usage()
        write(error_unit,'(A)') 'Usage: update_netgen_rate_from_dat DAT_FILE NETGEN_FILE REACTION [OUTPUT_FILE]'
        write(error_unit,'(A)') ''
        write(error_unit,'(A)') 'If OUTPUT_FILE is omitted, writes to tmp/<netgen-basename>.updated.txt.'
        write(error_unit,'(A)') 'If OUTPUT_FILE starts with /tmp/, it is redirected to project-local tmp/.'
        write(error_unit,'(A)') 'Grid points below the first finite article rate are log-log extrapolated downward.'
        write(error_unit,'(A)') ''
        write(error_unit,'(A)') 'Example:'
        write(error_unit,'(A)') '  update_netgen_rate_from_dat \'
        write(error_unit,'(A)') '    ppn_nova/references/iliadis2001_table06.dat \'
        write(error_unit,'(A)') '    NPDATA/netgen/netgen_iliadis2001_log_100.txt \'
        write(error_unit,'(A)') '    30P_pg_31S'
    end subroutine print_usage

    character(len=line_len) function normalized_output_path(requested_path, netgen_path) result(path)
        character(len=*), intent(in) :: requested_path, netgen_path
        character(len=line_len) :: name

        if (len_trim(requested_path) == 0) then
            name = path_basename(netgen_path)
            path = 'tmp/' // trim(updated_filename(name))
        else if (starts_with(trim(requested_path), '/tmp/')) then
            name = path_basename(requested_path)
            path = 'tmp/' // trim(name)
        else
            path = trim(requested_path)
        end if
    end function normalized_output_path

    character(len=line_len) function updated_filename(name) result(output)
        character(len=*), intent(in) :: name
        integer :: n

        n = len_trim(name)
        if (n > 4 .and. name(n-3:n) == '.txt') then
            output = name(1:n-4) // '.updated.txt'
        else
            output = trim(name) // '.updated.txt'
        end if
    end function updated_filename

    subroutine ensure_project_tmp_dir(path)
        character(len=*), intent(in) :: path
        integer :: status

        if (starts_with(trim(path), 'tmp/')) then
            call execute_command_line('mkdir -p tmp', exitstat=status)
            if (status /= 0) call fatal('could not create project tmp directory')
        end if
    end subroutine ensure_project_tmp_dir

    character(len=line_len) function path_basename(path) result(name)
        character(len=*), intent(in) :: path
        integer :: i

        name = trim(path)
        do i = len_trim(path), 1, -1
            if (path(i:i) == '/') then
                name = path(i + 1:len_trim(path))
                return
            end if
        end do
    end function path_basename

    subroutine read_article_rates(path, reaction_key, t9, rates, nsource)
        character(len=*), intent(in) :: path, reaction_key
        real(real64), allocatable, intent(out) :: t9(:), rates(:)
        integer, intent(out) :: nsource

        character(len=line_len) :: line
        character(len=token_len) :: tokens(max_tokens)
        integer :: unit, ios, ntokens, column_index, i
        real(real64) :: row_t9, row_rate
        logical :: have_header

        allocate(t9(max_points), rates(max_points))
        nsource = 0
        column_index = 0
        have_header = .false.

        open(newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) call fatal('could not open article dat file: ' // trim(path))

        do
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle

            call split_tokens(line, tokens, ntokens)
            if (ntokens == 0) cycle

            if (.not. have_header) then
                if (trim(tokens(1)) /= 'T9_GK') cycle
                do i = 2, ntokens
                    if (trim(canonical_key(tokens(i))) == trim(reaction_key)) then
                        column_index = i
                        exit
                    end if
                end do
                if (column_index == 0) then
                    call fatal('reaction column not found in dat file: ' // trim(reaction_key))
                end if
                have_header = .true.
                cycle
            end if

            if (ntokens < column_index) cycle
            if (trim(tokens(column_index)) == '...') cycle

            read(tokens(1),*,iostat=ios) row_t9
            if (ios /= 0) cycle
            read(tokens(column_index),*,iostat=ios) row_rate
            if (ios /= 0) cycle
            if (row_rate <= 0.0_real64) cycle

            nsource = nsource + 1
            if (nsource > max_points) call fatal('too many article rate points')
            t9(nsource) = row_t9
            rates(nsource) = row_rate
        end do

        close(unit)
        if (.not. have_header) call fatal('T9_GK header not found in dat file: ' // trim(path))
        if (nsource < 2) call fatal('need at least two finite positive source rates')
    end subroutine read_article_rates

    subroutine read_text_file(path, lines, nlines)
        character(len=*), intent(in) :: path
        character(len=line_len), allocatable, intent(out) :: lines(:)
        integer, intent(out) :: nlines

        character(len=line_len) :: line
        integer :: unit, ios

        nlines = 0
        open(newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) call fatal('could not open netgen file: ' // trim(path))
        do
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) exit
            nlines = nlines + 1
        end do
        close(unit)

        allocate(lines(nlines))
        open(newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) call fatal('could not reopen netgen file: ' // trim(path))
        do ios = 1, nlines
            read(unit,'(A)') lines(ios)
        end do
        close(unit)
    end subroutine read_text_file

    subroutine write_text_file(path, lines, nlines)
        character(len=*), intent(in) :: path
        character(len=line_len), intent(in) :: lines(:)
        integer, intent(in) :: nlines
        integer :: unit, ios, i

        open(newunit=unit, file=path, status='replace', action='write', iostat=ios)
        if (ios /= 0) call fatal('could not open output file: ' // trim(path))
        do i = 1, nlines
            write(unit,'(A)') trim(lines(i))
        end do
        close(unit)
    end subroutine write_text_file

    subroutine find_netgen_rate_block(lines, nlines, reaction_key, rate_start, rate_end)
        character(len=line_len), intent(in) :: lines(:)
        integer, intent(in) :: nlines
        character(len=*), intent(in) :: reaction_key
        integer, intent(out) :: rate_start, rate_end

        character(len=label_len) :: species_window(4), species
        character(len=label_len) :: block_key
        integer :: i, j
        logical :: matched_species

        species_window = ''
        rate_start = 0
        rate_end = 0

        do i = 1, nlines
            call parse_species_line(lines(i), matched_species, species)
            if (matched_species) then
                species_window(1:3) = species_window(2:4)
                species_window(4) = species
                cycle
            end if

            if (starts_with(trim(adjustl(lines(i))), '#Qrad')) then
                block_key = netgen_block_key(species_window)
                if (trim(block_key) == trim(reaction_key)) then
                    do j = i + 1, nlines
                        if (index(lines(j), '#       T8') > 0) exit
                    end do
                    if (j > nlines) call fatal('found reaction block but not T8 header')
                    rate_start = first_data_line(lines, nlines, j + 1)
                    if (rate_start == 0) call fatal('found reaction block but not rate table')
                    rate_end = rate_start
                    do while (rate_end <= nlines .and. is_rate_line(lines(rate_end)))
                        rate_end = rate_end + 1
                    end do
                    rate_end = rate_end - 1
                    return
                end if
            end if
        end do

        call fatal('reaction block not found in netgen file: ' // trim(reaction_key))
    end subroutine find_netgen_rate_block

    subroutine replace_rate_block(lines, rate_start, rate_end, source_t9, source_rates, nsource)
        character(len=line_len), intent(inout) :: lines(:)
        integer, intent(in) :: rate_start, rate_end, nsource
        real(real64), intent(in) :: source_t9(:), source_rates(:)

        real(real64), allocatable :: grid_t8(:)
        real(real64) :: t8, old_rate, t9, new_rate
        integer :: ngrid, i, ios

        ngrid = rate_end - rate_start + 1
        allocate(grid_t8(ngrid))

        do i = 1, ngrid
            read(lines(rate_start + i - 1),*,iostat=ios) t8, old_rate
            if (ios /= 0) call fatal('could not parse netgen rate line')
            grid_t8(i) = t8
        end do

        do i = 1, ngrid
            t8 = grid_t8(i)
            t9 = t8 / 10.0_real64
            if (t9 < source_t9(1)) then
                new_rate = log_extrapolate_first(source_t9, source_rates, nsource, t9)
            else
                new_rate = log_interpolate(source_t9, source_rates, nsource, t9)
            end if
            write(lines(rate_start + i - 1),'(F11.4,2X,A)') t8, trim(netgen_rate_string(new_rate))
        end do
    end subroutine replace_rate_block

    real(real64) function log_interpolate(t9_values, rate_values, nsource, t9) result(rate)
        real(real64), intent(in) :: t9_values(:), rate_values(:), t9
        integer, intent(in) :: nsource
        integer :: i
        real(real64) :: weight

        rate = -1.0_real64
        if (t9 < t9_values(1) .or. t9 > t9_values(nsource)) return

        do i = 1, nsource - 1
            if (t9 >= t9_values(i) .and. t9 <= t9_values(i + 1)) then
                if (abs(t9 - t9_values(i)) < 1.0e-14_real64) then
                    rate = rate_values(i)
                else if (abs(t9 - t9_values(i + 1)) < 1.0e-14_real64) then
                    rate = rate_values(i + 1)
                else
                    weight = (log10(t9) - log10(t9_values(i))) / &
                             (log10(t9_values(i + 1)) - log10(t9_values(i)))
                    rate = 10.0_real64 ** (log10(rate_values(i)) + &
                           weight * (log10(rate_values(i + 1)) - log10(rate_values(i))))
                end if
                return
            end if
        end do
    end function log_interpolate

    real(real64) function log_extrapolate_first(t9_values, rate_values, nsource, t9) result(rate)
        real(real64), intent(in) :: t9_values(:), rate_values(:), t9
        integer, intent(in) :: nsource
        real(real64) :: weight

        if (nsource < 2) then
            rate = -1.0_real64
            return
        end if

        weight = (log10(t9) - log10(t9_values(1))) / (log10(t9_values(2)) - log10(t9_values(1)))
        rate = 10.0_real64 ** (log10(rate_values(1)) + weight * (log10(rate_values(2)) - log10(rate_values(1))))
    end function log_extrapolate_first

    character(len=16) function netgen_rate_string(rate) result(out)
        real(real64), intent(in) :: rate
        real(real64) :: mantissa
        integer :: exponent
        character(len=8) :: mantissa_text
        character(len=4) :: exponent_text

        if (rate <= 0.0_real64) then
            out = '.9999E-99'
            return
        end if

        exponent = floor(log10(rate)) + 1
        mantissa = rate / (10.0_real64 ** exponent)
        if (mantissa >= 0.99995_real64) then
            mantissa = mantissa / 10.0_real64
            exponent = exponent + 1
        end if

        write(mantissa_text,'(F6.4)') mantissa
        if (exponent >= 0) then
            write(exponent_text,'(A,I2.2)') '+', exponent
        else
            write(exponent_text,'(A,I2.2)') '-', abs(exponent)
        end if
        out = mantissa_text(2:6) // 'E' // exponent_text
    end function netgen_rate_string

    integer function first_data_line(lines, nlines, start_index) result(index_out)
        character(len=line_len), intent(in) :: lines(:)
        integer, intent(in) :: nlines, start_index
        integer :: i

        index_out = 0
        do i = start_index, nlines
            if (is_rate_line(lines(i))) then
                index_out = i
                return
            end if
            if (starts_with(trim(adjustl(lines(i))), '#Qrad')) return
        end do
    end function first_data_line

    logical function is_rate_line(line) result(ok)
        character(len=*), intent(in) :: line
        character(len=line_len) :: clean
        real(real64) :: x, y
        integer :: ios

        clean = adjustl(line)
        if (len_trim(clean) == 0) then
            ok = .false.
            return
        end if
        if (clean(1:1) == '#') then
            ok = .false.
            return
        end if
        read(clean,*,iostat=ios) x, y
        ok = (ios == 0)
    end function is_rate_line

    subroutine parse_species_line(line, matched, species)
        character(len=*), intent(in) :: line
        logical, intent(out) :: matched
        character(len=label_len), intent(out) :: species

        character(len=line_len) :: clean
        character(len=token_len) :: symbol, mass
        integer :: count, ios

        matched = .false.
        species = ''

        clean = adjustl(line)
        if (len_trim(clean) == 0) return
        if (clean(1:1) /= '#') return

        clean = adjustl(clean(2:))
        symbol = ''
        mass = ''
        read(clean,*,iostat=ios) count, symbol, mass
        if (ios /= 0) then
            read(clean,*,iostat=ios) count, symbol
            if (ios /= 0) return
        end if

        matched = .true.
        symbol = upper_string(trim(symbol))

        if (count == 0 .or. trim(symbol) == 'OOOOO') then
            species = ''
        else if (trim(symbol) == 'PROT') then
            species = 'p'
        else
            species = isotope_label(trim(symbol), trim(mass))
        end if
    end subroutine parse_species_line

    character(len=label_len) function isotope_label(symbol, mass) result(label)
        character(len=*), intent(in) :: symbol, mass
        character(len=8) :: element
        character(len=16) :: digits, suffix
        integer :: i, nd, ns

        element = element_symbol(symbol)
        digits = ''
        suffix = ''
        nd = 0
        ns = 0

        do i = 1, len_trim(mass)
            if (mass(i:i) >= '0' .and. mass(i:i) <= '9') then
                nd = nd + 1
                digits(nd:nd) = mass(i:i)
            else
                ns = ns + 1
                suffix(ns:ns) = mass(i:i)
            end if
        end do

        label = trim(digits) // trim(element) // trim(suffix)
    end function isotope_label

    character(len=8) function element_symbol(symbol) result(element)
        character(len=*), intent(in) :: symbol

        select case (trim(upper_string(symbol)))
        case ('HE')
            element = 'He'
        case ('LI')
            element = 'Li'
        case ('BE')
            element = 'Be'
        case ('NE')
            element = 'Ne'
        case ('NA')
            element = 'Na'
        case ('MG')
            element = 'Mg'
        case ('AL')
            element = 'Al'
        case ('SI')
            element = 'Si'
        case ('CL')
            element = 'Cl'
        case ('AR')
            element = 'Ar'
        case ('CA')
            element = 'Ca'
        case ('SC')
            element = 'Sc'
        case default
            element = trim(symbol)
        end select
    end function element_symbol

    character(len=label_len) function netgen_block_key(species_window) result(key)
        character(len=label_len), intent(in) :: species_window(4)
        character(len=8) :: reaction_type

        key = ''
        if (len_trim(species_window(1)) == 0) return
        if (trim(species_window(2)) /= 'p') return
        if (len_trim(species_window(4)) == 0) return

        if (len_trim(species_window(3)) == 0) then
            reaction_type = 'pg'
        else if (trim(canonical_key(species_window(3))) == '4he') then
            reaction_type = 'pa'
        else
            return
        end if

        key = canonical_key(trim(species_window(1)) // '_' // trim(reaction_type) // '_' // trim(species_window(4)))
    end function netgen_block_key

    character(len=label_len) function canonical_key(raw) result(key)
        character(len=*), intent(in) :: raw

        character(len=line_len) :: lower
        integer :: i, n, position

        key = ''
        lower = lower_string(raw)
        i = 1
        n = len_trim(raw)
        position = 0

        do while (i <= n)
            if (i + 4 <= n .and. lower(i:i+4) == '(p,g)') then
                call append_text(key, position, '_pg_')
                i = i + 5
            else if (i + 4 <= n .and. lower(i:i+4) == '(p,a)') then
                call append_text(key, position, '_pa_')
                i = i + 5
            else if (i + 4 <= n .and. lower(i:i+4) == '(a,g)') then
                call append_text(key, position, '_ag_')
                i = i + 5
            else if (i + 4 <= n .and. lower(i:i+4) == '(a,p)') then
                call append_text(key, position, '_ap_')
                i = i + 5
            else if (raw(i:i) == ' ' .or. raw(i:i) == achar(9)) then
                i = i + 1
            else
                call append_text(key, position, lower(i:i))
                i = i + 1
            end if
        end do
    end function canonical_key

    subroutine split_tokens(line, tokens, ntokens)
        character(len=*), intent(in) :: line
        character(len=token_len), intent(out) :: tokens(:)
        integer, intent(out) :: ntokens

        integer :: i, n, start_pos
        logical :: in_token

        tokens = ''
        ntokens = 0
        n = len_trim(line)
        in_token = .false.
        start_pos = 1

        do i = 1, n + 1
            if (i <= n .and. line(i:i) /= ' ' .and. line(i:i) /= achar(9)) then
                if (.not. in_token) then
                    in_token = .true.
                    start_pos = i
                end if
            else
                if (in_token) then
                    ntokens = ntokens + 1
                    if (ntokens > size(tokens)) call fatal('too many tokens in line')
                    tokens(ntokens) = line(start_pos:i-1)
                    in_token = .false.
                end if
            end if
        end do
    end subroutine split_tokens

    subroutine append_text(buffer, position, text)
        character(len=*), intent(inout) :: buffer
        integer, intent(inout) :: position
        character(len=*), intent(in) :: text
        integer :: n

        n = len_trim(text)
        if (position + n > len(buffer)) call fatal('internal string buffer too short')
        buffer(position + 1:position + n) = text(1:n)
        position = position + n
    end subroutine append_text

    logical function starts_with(text, prefix) result(ok)
        character(len=*), intent(in) :: text, prefix
        if (len_trim(text) < len_trim(prefix)) then
            ok = .false.
        else
            ok = text(1:len_trim(prefix)) == prefix(1:len_trim(prefix))
        end if
    end function starts_with

    character(len=len(input)) function lower_string(input) result(output)
        character(len=*), intent(in) :: input
        integer :: i, c

        output = input
        do i = 1, len(input)
            c = iachar(input(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) output(i:i) = achar(c + 32)
        end do
    end function lower_string

    character(len=len(input)) function upper_string(input) result(output)
        character(len=*), intent(in) :: input
        integer :: i, c

        output = input
        do i = 1, len(input)
            c = iachar(input(i:i))
            if (c >= iachar('a') .and. c <= iachar('z')) output(i:i) = achar(c - 32)
        end do
    end function upper_string

    subroutine fatal(message)
        character(len=*), intent(in) :: message
        write(error_unit,'(A)') 'error: ' // trim(message)
        stop 1
    end subroutine fatal

end program update_netgen_rate_from_dat
