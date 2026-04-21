![](Resources/CAN_Satellite_2026_Poster.png)

# **Prerequisites and Installation**

To get started, you will need to install the required software and download the necessary data. It is important to install R before installing RStudio so that RStudio can automatically detect your R installation.

- Download R first from the CRAN website (version 4.2 or higher is recommended) [here](https://muug.ca/mirror/cran/).
Be sure to pick the R version specific to your operating system (Windows, macOS, or Linux),
and install the software in the default directories prompted by the installer.
- Download RStudio Desktop only after R is fully installed [here](https://posit.co/download/rstudio-desktop/https://posit.co/download/rstudio-desktop/).
Again, pick the RStudio version specific to your operating system,
and install it in the default directory.

- Ensure your system has at least 8GB of memory, though 16GB or more is highly recommended.

- Ensure atleast 15 GB of free storage space.

- Download the required UBC data objects (.rds files) from this OSF link under the Data folder [here](https://osf.io/rgajd/files/osfstorage).
Unzip these files in your downloads folder. These data were published and shared by M. Kapustina, details [here](https://osf.io/rgajd/files/6gxy3)

---

## **Overview**

Repository structure:
```
└── xenium_demo/
    ├── Readme.md  # Current file  
    ├── scripts/ 
    │   ├── Walkthrough.md
    │   ├── 0_Load_Packages.R
    |   ├── 1_Single_file_Demo.R
    │   ├── 2_Merging_Demo.R
    │   ├── 3_SingleR_lavel_transfer.R
    │   ├── 4_Differential_expression.R
    |
```


---
## **Open in Rstudio**

#### To download and open the repo in Rstudio, do the following: 

##### First download the repository:

- Click on the green button saying "code" at the top of this page
- Then Click on Download Zip 
- In your downloads folder, unzip the file in the downloads folder itself
(Note: Macs automatically unzip these files.)
 
##### Then Open Rstudio and:

- Click on file (top-left)
- Click on new project, then existing directory
- Navigate to your downloads folder and select the unzipped directory

##### Then, to begin demo:

- Navigate to the scripts folder
- Follow the scripts in numbered order (0 to 4)

> NOTE: Make sure to run 0_Load_Packages.R and download the data 
>       before the workshop begins.





