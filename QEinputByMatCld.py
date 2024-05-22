'''
conda env list
# conda environments:
#
base                  *  /opt/miniconda3
deepmd-cpu               /opt/miniconda3/envs/deepmd-cpu
pymatgen                 /opt/miniconda3/envs/pymatgen

conda remove --yes --name pymatgen --all
conda create --yes --name pymatgen python
source activate pymatgen  # OSX or Linux
conda install --yes --channel conda-forge pymatgen
conda install --yes numpy scipy matplotlib selenium
pip install --yes pymatgen selenium
'''

import time
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import NoSuchElementException, TimeoutException
import os


file_CWD = os.getcwd()

input_file_name = "input.in"
#input_file_name = "input.cif"

chrome_driver_path = "/opt/webdriver/chromedriver"

options = webdriver.ChromeOptions()
options.add_argument('--headless')
options.add_argument('--disable-gpu')
options.add_argument('--no-sandbox')

# 啟動Chrome WebDriver
driver = webdriver.Chrome(chrome_driver_path, options=options)

# 輸入檔案的路徑
input_file_path = os.path.join(os.getcwd(), input_file_name)

# 目標網站
driver.get("https://www.materialscloud.org/work/tools/qeinputgenerator")

#driver.implicitly_wait(5) #隱式等待，網頁載入分2階段，因為要等frame(0)出現而不是單純等待網頁載入所以不能用!!!
WebDriverWait(driver, 5).until(EC.frame_to_be_available_and_switch_to_it(0)) #顯式等待，把frame(0)出現並可切換當作等到條件
#------------------------------------
# 上傳結構文件
#driver.switch_to.frame(0)
Structurefile = driver.find_element(By.XPATH, "/html/body/div[1]/div[3]/form/div/div[1]/div[2]/input").send_keys(input_file_path)
driver.implicitly_wait(1)
#------------------------------------

#----------選FileformatSelect----------
fileformat = driver.find_element(By.ID, "fileformatSelect")
fileformat.click()
fileformat_Select = driver.find_element(By.XPATH, "/html/body/div[1]/div[3]/form/div/div[2]/div[2]/select/option[1]")### Quantum ESPRESSO input [parser: qetools]
#fileformat_Select = driver.find_element(By.XPATH, "/html/body/div[1]/div[3]/form/div/div[2]/div[2]/select/option[2]")### CIF File (.cif) [parser: ase]
#fileformat_Select = driver.find_element(By.XPATH, "/html/body/div[1]/div[3]/form/div/div[2]/div[2]/select/option[3]")### CIF File (.cif) [parser: pymatgen]/html/body/div[1]/div[3]/form/div/div[2]/div[2]/select/option[3]
fileformat_Select.click()
#------------------------------------

#----------選PseudoSelect----------
Pseudo = driver.find_element(By.ID, "pseudoSelect")
Pseudo.click()
#time.sleep(1)
#Pseudo_Select = driver.find_element(By.XPATH, "/html/body/div[1]/div[3]/form/div/div[4]/div[2]/select/option[1]")### SSSP Efficiency PBEsol (version 1.1)
#Pseudo_Select = driver.find_element(By.XPATH, "/html/body/div[1]/div[3]/form/div/div[4]/div[2]/select/option[2]")### SSSP Precision PBEsol (version 1.1)
Pseudo_Select = driver.find_element(By.XPATH, "/html/body/div[1]/div[3]/form/div/div[4]/div[2]/select/option[3]")### SSSP Efficiency PBE (version 1.1)
#Pseudo_Select = driver.find_element(By.XPATH, "/html/body/div[1]/div[3]/form/div/div[4]/div[2]/select/option[4]")### SSSP Precision PBE (version 1.1)
Pseudo_Select.click()
#------------------------------------

#----------選MagnetizationSelect----------
Magnetization = driver.find_element(By.ID, "magnetizationSelect")
Magnetization.click()
#time.sleep(1)
#MagnetizationSelect = driver.find_element(By.XPATH, "/html/body/div[1]/div[3]/form/div/div[5]/div[2]/select/option[1]")### non-magnetic metal (fractional occupations)
#MagnetizationSelect = driver.find_element(By.XPATH, "/html/body/div[1]/div[3]/form/div/div[5]/div[2]/select/option[2]")### non-magnetic insulator (fixed occupations)
MagnetizationSelect = driver.find_element(By.XPATH, "/html/body/div[1]/div[3]/form/div/div[5]/div[2]/select/option[3]")### magnetic (fractional occupations)
MagnetizationSelect.click()
#------------------------------------

#---------選KmeshSelect-----------
Kmesh = driver.find_element(By.ID, "kmeshSelect")
Kmesh.click()
#time.sleep(1)
#KmeshSelect = driver.find_element(By.XPATH, "/html/body/div[1]/div[3]/form/div/div[6]/div[2]/select/option[1]")### very fine (0.15 1/Å, 0.1 eV)
#KmeshSelect = driver.find_element(By.XPATH, "/html/body/div[1]/div[3]/form/div/div[6]/div[2]/select/option[2]")### fine (0.20 1/Å, 0.2 eV)
KmeshSelect = driver.find_element(By.XPATH, "/html/body/div[1]/div[3]/form/div/div[6]/div[2]/select/option[3]")### normal (0.30 1/Å, 0.3 eV)
KmeshSelect.click()
#------------------------------------

#----------選RefineCellSelect----------
#RefineCell = driver.find_element(By.ID, "refineCellSelect")
#RefineCell.click()
#time.sleep(1)
#RefineCellSelect = driver.find_element(By.XPATH, "/html/body/div[1]/div[3]/form/div/div[7]/div[2]/select/option[1]")### No
## RefineCellSelect = driver.find_element(By.XPATH, "/html/body/div[1]/div[3]/form/div/div[7]/div[2]/select/option[2]")### Yes (choose the symmetry tolerance below)
#RefineCellSelect.click()
#------------------------------------

#------------------------------------
###----------按Generate the PWscf input file----------
Generate_the_PWscf_input_file = driver.find_element(By.XPATH, "/html/body/div[1]/div[3]/form/div/div[8]/input")
Generate_the_PWscf_input_file.click()
#------------------------------------
WebDriverWait(driver, 5)
#------------------------------------
##------有Warning按Got it!------------------
try:
    button = WebDriverWait(driver, 2).until(
            EC.presence_of_element_located((By.XPATH, "/html/body/div[4]/div/div/div[3]/button"))
    )
    do_generate_pwscf_input_file = True
except (NoSuchElementException,TimeoutException):
    do_generate_pwscf_input_file = False
if do_generate_pwscf_input_file:
    Generate_the_PWscf_input_file = driver.find_element(By.XPATH, "/html/body/div[4]/div/div/div[3]/button")
    Generate_the_PWscf_input_file.click()
#    driver.implicitly_wait(5)
    
# 拿PWscf_input
PWscf_inputs = driver.find_element(By.ID, "qepwinput")
text = PWscf_inputs.text

# 寫入output.in
output_file_name = "output.in"
with open(output_file_name, 'w') as output_file:
    output_file.write(text)

# 關掉chrome
driver.quit()
#print("PWscf get!")
