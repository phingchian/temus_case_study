FROM mambaorg/micromamba:1.5.8

LABEL image.author.name="Chai Phing Chian"
LABEL image.author.email="imcpc@pm.me"

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml .

RUN micromamba install -y -f env.yaml && \
    micromamba clean --all --yes

ENV PATH=/opt/conda/bin:/tmp:$PATH

RUN python -m wget "https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20240625.zip"

RUN unzip plink2_linux_x86_64_20240625.zip \
    && rm plink2_linux_x86_64_20240625.zip