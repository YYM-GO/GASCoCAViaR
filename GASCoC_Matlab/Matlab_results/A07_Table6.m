close
clear
clc

load shenwan
load hushen
load covarportGauSAV
load coesportGauSAV
load weightGauSAV
load covarportGauAS
load coesportGauAS
load weightGauAS
load covarportPOTSAV
load coesportPOTSAV
load weightPOTSAV
load covarportPOTAS
load coesportPOTAS
load weightPOTAS

load covarportGauSAVc
load coesportGauSAVc
load weightGauSAVc
load covarportGauASc
load coesportGauASc
load weightGauASc
load covarportPOTSAVc
load coesportPOTSAVc
load weightPOTSAVc
load covarportPOTASc
load coesportPOTASc
load weightPOTASc

load weightcopulaGauall
load copulaCoESGauport
load copulaCoVaRGauport

load weightcopulaPOTall
load copulaCoVaRPOTport
load copulaCoESPOTport

load weightDCCGauall
load DCCCoVaRGauport
load DCCCoESGauport

load weightDCCPOTall
load DCCCoVaRPOTport
load DCCCoESPOTport

[T,N] = size(shenwan);
wind = 1400;
Tout = T-wind;

%%                        Panel A
portfolio_GARCH_Gaussian=[mean(covarportGauSAV),std(covarportGauSAV),mean(coesportGauSAV),std(coesportGauSAV),skewport(weightGauSAV,shenwan(wind+1:end,:)),kurport(weightGauSAV,shenwan(wind+1:end,:));
    mean(covarportGauAS),std(covarportGauAS),mean(coesportGauAS),std(coesportGauAS),skewport(weightGauAS,shenwan(wind+1:end,:)),kurport(weightGauAS,shenwan(wind+1:end,:));
    mean(covarportGauSAVc),std(covarportGauSAVc),mean(coesportGauSAVc),std(coesportGauSAVc),skewport(weightGauSAVc,shenwan(wind+1:end,:)),kurport(weightGauSAVc,shenwan(wind+1:end,:));
    mean(covarportGauASc),std(covarportGauASc),mean(coesportGauASc),std(coesportGauASc),skewport(weightGauASc,shenwan(wind+1:end,:)),kurport(weightGauASc,shenwan(wind+1:end,:));
    mean(copulaCoVaRGauport(:,3)),std(copulaCoVaRGauport(:,3)),mean(copulaCoESGauport(:,3)),std(copulaCoESGauport(:,3)),skewport(weightcopulaGauall,shenwan(wind+1:end,:)),kurport(weightcopulaGauall,shenwan(wind+1:end,:));
    mean(copulaCoVaRGauport(:,2)),std(copulaCoVaRGauport(:,2)),mean(copulaCoESGauport(:,2)),std(copulaCoESGauport(:,2)),skewport(weightcopulaGauall,shenwan(wind+1:end,:)),kurport(weightcopulaGauall,shenwan(wind+1:end,:));
    mean(copulaCoVaRGauport(:,4)),std(copulaCoVaRGauport(:,4)),mean(copulaCoESGauport(:,4)),std(copulaCoESGauport(:,4)),skewport(weightcopulaGauall,shenwan(wind+1:end,:)),kurport(weightcopulaGauall,shenwan(wind+1:end,:));
    mean(copulaCoVaRGauport(:,1)),std(copulaCoVaRGauport(:,1)),mean(copulaCoESGauport(:,1)),std(copulaCoESGauport(:,1)),skewport(weightcopulaGauall,shenwan(wind+1:end,:)),kurport(weightcopulaGauall,shenwan(wind+1:end,:));
    mean(DCCCoVaRGauport),std(DCCCoVaRGauport),mean(DCCCoESGauport),std(DCCCoESGauport),skewport(weightDCCGauall,shenwan(wind+1:end,:)),kurport(weightDCCGauall,shenwan(wind+1:end,:))]


%%                        Panel B
portfolio_GARCH_POT=[mean(covarportPOTSAV),std(covarportPOTSAV),mean(coesportPOTSAV),std(coesportPOTSAV),skewport(weightPOTSAV,shenwan(wind+1:end,:)),kurport(weightPOTSAV,shenwan(wind+1:end,:));
    mean(covarportPOTAS),std(covarportPOTAS),mean(coesportPOTAS),std(coesportPOTAS),skewport(weightPOTAS,shenwan(wind+1:end,:)),kurport(weightPOTAS,shenwan(wind+1:end,:));
    mean(covarportPOTSAVc),std(covarportPOTSAVc),mean(coesportPOTSAVc),std(coesportPOTSAVc),skewport(weightPOTSAVc,shenwan(wind+1:end,:)),kurport(weightPOTSAVc,shenwan(wind+1:end,:));
    mean(covarportPOTASc),std(covarportPOTASc),mean(coesportPOTASc),std(coesportPOTASc),skewport(weightPOTASc,shenwan(wind+1:end,:)),kurport(weightPOTASc,shenwan(wind+1:end,:));
    mean(copulaCoVaRPOTport(:,3)),std(copulaCoVaRPOTport(:,3)),mean(copulaCoESPOTport(:,3)),std(copulaCoESPOTport(:,3)),skewport(weightcopulaPOTall,shenwan(wind+1:end,:)),kurport(weightcopulaPOTall,shenwan(wind+1:end,:));
    mean(copulaCoVaRPOTport(:,2)),std(copulaCoVaRPOTport(:,2)),mean(copulaCoESPOTport(:,2)),std(copulaCoESPOTport(:,2)),skewport(weightcopulaPOTall,shenwan(wind+1:end,:)),kurport(weightcopulaPOTall,shenwan(wind+1:end,:));
    mean(copulaCoVaRPOTport(:,4)),std(copulaCoVaRPOTport(:,4)),mean(copulaCoESPOTport(:,4)),std(copulaCoESPOTport(:,4)),skewport(weightcopulaPOTall,shenwan(wind+1:end,:)),kurport(weightcopulaPOTall,shenwan(wind+1:end,:));
    mean(copulaCoVaRPOTport(:,1)),std(copulaCoVaRPOTport(:,1)),mean(copulaCoESPOTport(:,1)),std(copulaCoESPOTport(:,1)),skewport(weightcopulaPOTall,shenwan(wind+1:end,:)),kurport(weightcopulaPOTall,shenwan(wind+1:end,:));
    mean(DCCCoVaRPOTport),std(DCCCoVaRPOTport),mean(DCCCoESPOTport),std(DCCCoESPOTport),skewport(weightDCCPOTall,shenwan(wind+1:end,:)),kurport(weightDCCPOTall,shenwan(wind+1:end,:))]

