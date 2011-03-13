#Set release name or install directory
RELEASE_NAME=mugsy_annotator-0.5
INSTALL_DIR=./${RELEASE_NAME}

all: mugsya_install

install: mugsya_install

dist:
	tar cvzf ${RELEASE_NAME}.tgz ${INSTALL_DIR}

mugsya_install:
	mkdir -p ${INSTALL_DIR}
	install README ${INSTALL_DIR}
	install mugsy-annotator ${INSTALL_DIR}
	install mapfeatures.pl ${INSTALL_DIR}
	install mafindex.pl ${INSTALL_DIR}
	install featureindex.pl ${INSTALL_DIR}
	install AlignmentTree.pm ${INSTALL_DIR}
	install IntervalTree.pm ${INSTALL_DIR}
