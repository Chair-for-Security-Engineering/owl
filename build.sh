#!/bin/sh

CFLAGS="-O3 -march=native -Wall -Wconversion -Wextra -Wpedantic -Wshadow -Wundef -Wunused"
CXX="g++ -std=c++17 ${CFLAGS} -fpermissive"
CC="gcc -std=c99 ${CFLAGS}"

SSHOPT="-o BatchMode=yes -o ConnectTimeout=2"
SCP="scp ${SSHOPT} -r"
SSH="ssh ${SSHOPT}"


build() {
	echo "\033[36m./build.sh: build\033[0m"
	incs="-I'${2}/hexl/hexl/include'"
	libs="-L${2}/build/hexl/hexl/lib -lhexl"
	${SSH} "${1}" "{
		if ! test -d '${2}/build/hexl'; then
			cmake -S '${2}/hexl' -B '${2}/build/hexl' && \
			cmake --build '${2}/build/hexl' -j
		fi && \
		${CC}  -c -o '${2}/build/test.o' '${2}/test.c' && \
		${CXX} -c -o '${2}/build/ring.o' '${2}/ring_hexl.cpp' ${incs} && \
		${CXX}    -o '${2}/build/test'   '${2}/build/'{test,ring}.o ${libs}
	}" || exit 1
}

clean() {
	echo "\033[36m./build.sh: clean\033[0m"
	${SSH} "${1}" "rm -r '${2}'" || exit 1
}

rcopy() {
	echo "\033[36m./build.sh: rcopy\033[0m"
	${SSH} "${1}" "mkdir -p '${2}'" || exit 1

	${SSH} "${1}" "test -d '${2}/hexl'" || {
		${SCP} hexl/* "${1}:${2}/hexl/" > /dev/null || exit 1
	}
	${SCP} build.sh ring.h ring_hexl.cpp swk.c test.c "${1}:${2}/" > /dev/null || exit 1
}

tests() {
	echo "\033[36m./build.sh: tests\033[0m"
	${SSH} "${1}" "{
		${2}/build/test
	}" || exit 1
}

usage() {
	echo "./build.sh [-h] [-r remote] [-d dir] command"  2>&1
	echo "  -h: print help and exit"                     2>&1
	echo "  -r: remote server (default: \$BUILD_REMOTE)" 2>&1
	echo "  -d: remote directory (default: \$BUILD_DIR)" 2>&1
	echo "\ncommands:"                                   2>&1
	echo "  build: build on remote"                      2>&1
	echo "  clean: remove build directory on remote"     2>&1
	echo "  rcopy: copy files to remote"                 2>&1
	echo "  tests: run tests on remote"                  2>&1
}


dir="${BUILD_DIR}"
remote="${BUILD_REMOTE}"
while getopts "d:hr:" opt; do
	case "${opt}" in
	d) dir="${OPTARG}";;
	h) usage; exit 0;;
	r) remote="${OPTARG}";;
	*) usage; exit 1;;
	esac
done
shift $(( OPTIND - 1 ))

test ${#} -eq 0 && { usage; exit 1; }
while test ${#} -gt 0; do
	case "${1}" in
	build) shift; build "${remote}" "${dir}";;
	clean) shift; clean "${remote}" "${dir}";;
	rcopy) shift; rcopy "${remote}" "${dir}";;
	tests) shift; tests "${remote}" "${dir}";;
	*)     usage; exit 1;
	esac
done
