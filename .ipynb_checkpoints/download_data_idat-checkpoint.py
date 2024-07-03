import os
import pathlib

gdc_client = '/media/hieunguyen/HNSD01/gdc-client'
infopath = "/media/hieunguyen/HNSD01/src/TCGA_data_analysis/TCGA_database_DMR"
savedir = "/media/hieunguyen/HNHD01/TCGA_all_idat/20240619"
tracking_file = "/media/hieunguyen/HNSD01/src/TCGA_data_analysis/finished_download_files.txt"

all_idat_folders = [item for item in pathlib.Path(infopath).glob("*idat*")]

for folder in all_idat_folders:
    label = folder.name.replace("_idat", "")
    os.system("mkdir -p {}".format(os.path.join(savedir, label, "tumor")))
    os.system("mkdir -p {}".format(os.path.join(savedir, label, "normal")))
    
    normal_manifest = [item for item in pathlib.Path(folder).glob("*manifest*Normal*450K*")]
    tumor_manifest = [item for item in pathlib.Path(folder).glob("*manifest*Tumor*450K*")]
    
    if len(normal_manifest) == 1:
        normal_manifest = normal_manifest[0]
        downloaddir = os.path.join(savedir, "normal")
        os.system("{} download --manifest {} --dir {}".format(gdc_client, str(normal_manifest), downloaddir))
        os.system("echo {} normal >> {}".format(label, tracking_file))
    if len(tumor_manifest) == 1:
        tumor_manifest = tumor_manifest[0]
        downloaddir = os.path.join(savedir, "tumor")
        os.system("{} download --manifest {} --dir {}".format(gdc_client, str(tumor_manifest), downloaddir))
        os.system("echo {} tumor >> {}".format(label, tracking_file))