#Set release name or install directory
RELEASE_NAME=mugsy_x86-64-v1r2.3.1
INSTALL_DIR=./${RELEASE_NAME}


all: nucmer synchain_mugsy mugsy_seqan

install: mugsy_install mummer_install

dist:
	tar cvzf ${RELEASE_NAME}.tgz ${INSTALL_DIR}
nucmer:
	make -C MUMmer3.20 all

synchain_mugsy:
	make -C chaining synchain-mugsy

mugsy_seqan:
	make -C mugsy-seqan Project=mugsy

mugsy_install:
	mkdir -p ${INSTALL_DIR}
	install mugsyenv.sh ${INSTALL_DIR}
	perl -pi -e 's|export MUGSY_INSTALL=.*|export MUGSY_INSTALL=${INSTALL_DIR}|' ${INSTALL_DIR}/mugsyenv.sh


	install mugsy ${INSTALL_DIR}
	install mugsyWGA ${INSTALL_DIR}
	install synchain-mugsy ${INSTALL_DIR}

	install maf2fasta.pl ${INSTALL_DIR}
	install xmfa2maf.pl ${INSTALL_DIR}
	install fixMAFnames.pl ${INSTALL_DIR}
	install splitmaf.pl ${INSTALL_DIR}
	install plot.pl ${INSTALL_DIR}
	install delta-dups.sh ${INSTALL_DIR}

	install -m 644 README ${INSTALL_DIR}
	install -m 644 CHANGELOG ${INSTALL_DIR}
	install -m 644 LICENSE ${INSTALL_DIR}
	install -m 644 LICENSE ${INSTALL_DIR}
	install -m 644 README.install ${INSTALL_DIR}

mummer_install:
	mkdir -p ${INSTALL_DIR}/MUMmer3.20/scripts

	install MUMmer3.20/delta-filter ${INSTALL_DIR}/MUMmer3.20/
	install MUMmer3.20/gaps ${INSTALL_DIR}/MUMmer3.20/
	install MUMmer3.20/mgaps ${INSTALL_DIR}/MUMmer3.20/
	install MUMmer3.20/delta2maf ${INSTALL_DIR}/MUMmer3.20/
	install MUMmer3.20/aux_bin/postnuc ${INSTALL_DIR}/MUMmer3.20/
	install MUMmer3.20/aux_bin/prenuc ${INSTALL_DIR}/MUMmer3.20/
#	install MUMmer3.20/src/tigr/show-coords ${INSTALL_DIR}/MUMmer3.20/
	install MUMmer3.20/mummer ${INSTALL_DIR}/MUMmer3.20/
	install MUMmer3.20/mummerplot ${INSTALL_DIR}/MUMmer3.20/
	install MUMmer3.20/nucmer ${INSTALL_DIR}/MUMmer3.20/
	install MUMmer3.20/scripts/Foundation.pm ${INSTALL_DIR}/MUMmer3.20/scripts

# mapping_install:
#	@install -d perllibs ${INSTALL_DIR}
#	@install mapping/mafindex.pl ${INSTALL_DIR}
#	@install mapping/AlignmentTree.pm ${INSTALL_DIR}
#	@install mapping/IntervalTree.pm ${INSTALL_DIR}
#	@install mapping/featureindex.pl ${INSTALL_DIR}
#	@install mapping/mapfeatures.pl ${INSTALL_DIR}

# multiz_install:
#	@install labelblocks.pl ${INSTALL_DIR}
