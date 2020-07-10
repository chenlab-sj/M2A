FROM python:3.6.5
RUN pip install --upgrade pip

WORKDIR /M2A

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY 1_ResponseVariable 1_ResponseVariable/
COPY 2_MethylationFeatures 2_MethylationFeatures/
COPY 3_CombineInput 3_CombineInput/
COPY 4_TransferLearning 4_TransferLearning/
COPY 5_RunModel 5_RunModel/

ENTRYPOINT ["/bin/bash"]
