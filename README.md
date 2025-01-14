Code and data associated with the manuscript: 

    Decoding Prosodic Information from Motion Capture Data: The Gravity of Co-Speech Gestures 

Citation: 

    @inproceedings{--,
      title={Decoding Prosodic Information from Motion Capture Data: The Gravity of Co-Speech Gestures},
      author={Momsen, Jacob P. and Coulson, Seana},
      journaltitle = {Open Mind},
      year={2025},
      month = {--},
      url = {https://},
      month_numeric = {--}
    }

For questions and issues, email jacob.momsen@yale.edu


Data used with permission from the Trinity Speech-Gesture Database (Database I) (https://trinityspeechgesture.scss.tcd.ie)
    Citation:
    
    @inproceedings{IVA:2018,
      title={IVA: Investigating the use of recurrent motion modelling for speech gesture generation},
      author={Ferstl, Ylva and McDonnell, Rachel},
      booktitle = {IVA '18 Proceedings of the 18th International Conference on Intelligent Virtual Agents},
      year={2018},
      month = {Nov},
      url = {https://trinityspeechgesture.scss.tcd.ie},
      month_numeric = {11}
    }

Primary script: 

    gocsg_main            
    
Functions:

    gocsg_TRFloop         # Runs the decoder analysis 
    gocsg_bvh2skel        # Transforms raw mocap data (.bvh) to MATLAB-compatible format

Folder stucture:
    
     ── Main                            # Primary folder - set as "homefolder" in gocsg_main
        ├── AVdata/                     
        │   ├── speech/                  
                │   ├── wav/            # Raw audio data from Trinity Speech-Gesture Database I (23 files)
                │   ├── envfiles/       # Stores speech data configured for TRF analysis 
        │   ├── mocap/                  
                │   ├── bvh/            # Raw mocap data from Trinity Speech-Gesture Database I (23 files)
                │   ├── skelfiles/      # Stores mocap data configured for TRF analysis 
        ├── setup_files/                # Contains required files for functions
        ├── TRFoutput/
        │   ├── Study1/                 # Stores output associated with Study 1 analysis 
        │   ├── Study2/                 # Stores output associated with Study 2 analysis
        ├── README.md                   # Overview of the repository
        ├── LICENSE                     # Licensing information

License:

    This project is licensed under the MIT License - see the LICENSE file for details.
    
Required external toolboxes (3): 
    
    https://github.com/mocaptoolbox/mocaptoolbox/tree/main
    https://github.com/mickcrosse/mTRF-Toolbox/tree/master
    https://github.com/alexisdmacintyre/AcousticLandmarks

