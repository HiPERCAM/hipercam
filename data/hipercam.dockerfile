FROM python:3.11-bookworm

# Installing pgplot5 and some other dependencies
RUN echo "deb http://deb.debian.org/debian bookworm main non-free" >> /etc/apt/sources.list
# Libqt5gui and gl1 are needed by Qt5Agg
RUN apt update && apt install -y pgplot5 gfortran libqt5gui5 libgl1-mesa-glx
# QT5 and dependencies are needed by PyQt5
#RUN apt install -y pyqt5-dev pyqt5-dev-tools libqt5core5a libqt5gui5 libqt5svg5 libqt5widgets5 python3-pyqt5
#ENV PYTHONPATH=/usr/lib/python3/dist-packages

# This step is required so that the cpgplot library is found by pgplot
RUN ln -s /usr/include/cpgplot.h /usr/lib/

# Upgrade pip and install some basic dependencies
RUN pip install --upgrade pip
RUN pip install Cython ipython
# pin numpy to >=2.0 because of changes to API
RUN pip install numpy>=2

# Installing PyQt5 
RUN pip install PyQt5
ENV MPLBACKEND=Qt5Agg

# Setting PGPLOT env variables
ENV PGPLOT_DIR=/usr/lib
ENV PGPLOT_PNG=true

# fixes libgl error with matplotlib
ENV LIBGL_ALWAYS_INDIRECT=1

# Working on the hiperuser directory but still as root, 
# so that pip dependencies are installed at container root level 
# and available alongside the preinstalled version of python
WORKDIR /home/hiperuser/software

# Pulling both required repositories
RUN git clone https://github.com/HiPERCAM/trm-pgplot.git && git clone https://github.com/HiPERCAM/hipercam.git

# Installing Tom's pgplot
WORKDIR /home/hiperuser/software/trm-pgplot
RUN pip install .

# Installing hipercam dependencies
WORKDIR /home/hiperuser/software/hipercam
RUN pip install .

# Installing hcam_finders and dependencies
WORKDIR /home/hiperuser
RUN pip install hcam_finder

# Now we add a non-priviledged user and give it rights in the /home/hiperuser folder.
RUN useradd -m hiperuser
RUN chown -R hiperuser:hiperuser /home/hiperuser
USER hiperuser

# Going back to the home user's directory. 
WORKDIR /home/hiperuser

CMD ["/bin/bash"]