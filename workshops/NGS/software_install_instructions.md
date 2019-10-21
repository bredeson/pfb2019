Install instructions for the software required by the NGS Workshop Problemset
==============================================

get the software precompiled via https://www.dropbox.com/s/zn4m0oqbasg4gal/software.tar.gz


FastQC
------
1. Download the [Mac DMG image](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.dmg) of FastQC.

2. Double-click the fastqc_v0.11.8.dmg to initialize the installation prompt.

3. A Finder window will appear with the FastQC logo. Click and drag this logo to the /Applications folder

4. In your favorite text editor, open `/Applications/FastQC.app/Contents/MacOS/fastqc` and change
   ```perl
   #!/usr/bin/perl
   ```
   To
   ```perl
   #!/usr/bin/env perl
   ```

5. Next, `cd /usr/local/bin` and `ln -s /Applications/FastQC.app/Contents/MacOS/fastqc .`, then `cd $HOME`.

6. Type `fastqc --help` to check that the software works.


SRA Toolkit
-----------
1. Download the SRA Toolkit [tarball](https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz)

2. Unpack the tarball `tar -xvf sratoolkit.current-mac64.tar`

3. Move the directory into place and link the executables into /usr/local/bin:
   ```bash
   mv sratoolkit.2.9.6-1-mac64 /usr/local
   pushd /usr/local/bin
   ln -s /usr/local/sratoolkit.2.9.6-1-mac64/bin/* .
   popd
   ```


