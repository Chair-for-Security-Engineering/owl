#!/bin/sh

export MY_CFLAGS="-O3 -march=native -Wall -Wconversion -Wextra -Wpedantic -Wshadow -Wundef -Wunused"
export MY_CXX="${CXX:-g++} -std=c++17 ${MY_CFLAGS} -fpermissive"
export MY_CC="${CC:-gcc} -std=c99 ${MY_CFLAGS}"

export MY_SSHOPT="-o BatchMode=yes -o ConnectTimeout=2"
export MY_SCP="scp ${MY_SSHOPT} -r"
export MY_SSH="ssh ${MY_SSHOPT}"

build() {
	if test -n "${1}"; then
		echo "\033[36m./build.sh: build\033[0m"
		${MY_SSH} "${1}" "${2}/build.sh build || exit 1" && return 0 || exit 1
	fi

	d="$(pwd)/$(dirname "${0}")"
	echo "\033[34m./build.sh: build\033[0m"

	b="${d}/build"
	mkdir -p "${b}" || exit 1

	${MY_CC} -c -o "${b}/test.o"  "${d}/test.c"  || exit 1
	${MY_CC} -c -o "${b}/bench.o" "${d}/bench.c" || exit 1

	if test "$(uname -m)" = "x86_64"; then
		if ! test -d "${b}/hexl"; then
			cmake -S "${d}/hexl" -B "${b}/hexl" || exit 1
			cmake --build "${b}/hexl" -j        || exit 1
		fi
		IHEXL="-I${d}/hexl/hexl/include"
		LHEXL="-L${b}/hexl/hexl/lib -lhexl -lm"
		${MY_CXX} -c -o "${b}/ring_hexl.o" "${d}/ring_hexl.cpp" ${IHEXL}              || exit 1
		${MY_CXX}    -o "${b}/test_hexl"   "${b}/ring_hexl.o" "${b}/test.o" ${LHEXL}  || exit 1
		${MY_CXX}    -o "${b}/bench_hexl"  "${b}/ring_hexl.o" "${b}/bench.o" ${LHEXL} || exit 1
	fi

	if test -f "${d}/latex.c"; then
		${MY_CC} -o "${b}/latex" "${d}/latex.c" -lm || exit 1
	fi

	return 0
}

clean() {
	if test -n "${1}"; then
		echo "\033[36m./build.sh: clean\033[0m"
		${MY_SSH} "${1}" "rm -rf ${2}" || exit 1
	else
		echo "\033[34m./build.sh: clean\033[0m"
		rm -rf "build" || exit 1
	fi

	return 0
}

rcopy() {
	if test -n "${1}"; then
		echo "\033[36m./build.sh: rcopy\033[0m"
		${MY_SSH} "${1}" "mkdir -p ${2}" || exit 1
		${MY_SSH} "${1}" "test -d ${2}/hexl" || {
			${MY_SCP} hexl/* "${1}:${2}/hexl/" > /dev/null || exit 1
		}
		FILES="build.sh ring.h ring_hexl.cpp swk.c test.c params.c bench.c latex.c"
		${MY_SCP} ${FILES} "${1}:${2}/" > /dev/null || exit 1
	fi

	return 0
}

tests() {
	if test -n "${1}"; then
		echo "\033[36m./build.sh: tests\033[0m"
		${MY_SSH} "${1}" "${2}/build.sh tests || exit 1" && return 0 || exit 1
	fi

	d="$(pwd)/$(dirname "${0}")"
	echo "\033[34m./build.sh: tests\033[0m"
	test -x "${d}/build/test_hexl" && { "${d}/build/test_hexl" || exit 1; }
	return 0
}

usage() {
	echo "build.sh [-h] [-r remote] command [command [...]]" 2>&1
	echo "  -h: print help and exit"                         2>&1
	echo "  -r: remote server (default: \"\")"               2>&1
	echo "\ncommands:"                                       2>&1
	echo "  rcopy: copy files to remote"                     2>&1
	echo "  build: build (on remote)"                        2>&1
	echo "  tests: run tests (on remote)"                    2>&1
	echo "  clean: remove directory (on remote)"             2>&1
}


dir="owl"
while getopts "hr:" opt; do
	case "${opt}" in
	h) usage; exit 0;;
	r) remote="${OPTARG}";;
	*) usage; exit 1;;
	esac
done
shift $(( OPTIND - 1 ))

test ${#} -eq 0 && { usage; exit 1; }
while test ${#} -gt 0; do
	case "${1}" in
	build) shift; build "${remote}" "${dir}" || exit 1;;
	clean) shift; clean "${remote}" "${dir}" || exit 1;;
	rcopy) shift; rcopy "${remote}" "${dir}" || exit 1;;
	tests) shift; tests "${remote}" "${dir}" || exit 1;;
	*)     usage; exit 1;
	esac
done
