FROM gcc
# Instructions copied from https://thenewstack.io/ruby-in-containers/
RUN curl -sSL https://rvm.io/mpapis.asc | gpg --import -
RUN curl -sSL https://rvm.io/pkuczynski.asc | gpg --import -
RUN curl -L https://get.rvm.io | bash -s stable
ENV PATH /usr/local/rvm/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
RUN rvm install ruby-2.7.1
COPY . /Lilly-Medchem-Rules
WORKDIR /Lilly-Medchem-Rules
RUN make
RUN [ "/bin/bash", "-l", "-c", ". /etc/profile && rvm use 2.7.1 > /dev/null && make test" ]
# docker run -i -a stdout -a stdin ianwatson/lilly_medchem_rules:v1.2 < file.smi
CMD /bin/bash -c ". /etc/profile && rvm use 2.7.1 > /dev/null && /Lilly-Medchem-Rules/Lilly_Medchem_Rules.rb -i smi -"
