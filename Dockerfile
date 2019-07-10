FROM kbase/sdkbase2:python
MAINTAINER Sean Jungbluth <sjungbluth@lbl.gov>
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.


# To install all the dependencies
RUN apt-get update && apt-get install -y libgsl0-dev git zip unzip python-pip && \
    apt-get install -y r-base r-cran-ggplot2

#RUN git clone https://github.com/jungbluth/kb_iMAG-viz

ADD https://github.com/jungbluth/kb_iMAG-viz/blob/master/make_ggplot_scatterplot.R /usr/local/bin/

RUN echo 'install.packages("ggpubr", dependencies=TRUE, repos = "http://cran.us.r-project.org")' > /tmp/packages.R \
    && Rscript /tmp/packages.R

#ENTRYPOINT [ "./scripts/entrypoint.sh" ]

#  write.csv(taxa, 'taxa.csv', col.names=NA); \n \
# RUN printf "library(ggplot2); \n \
#   library(ggpubr); \n \
#   \n" \
#   > ~/commands.R


WORKDIR /tmp
USER root


# # set the entrypoint
#ENTRYPOINT ["sh","-c","Rscript ~/commands.R"]

#CMD [ ]

