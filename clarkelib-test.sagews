︠2911106f-0cde-4040-91a2-48e8ba4df5b8is︠
%md
<font size=5> Load ClarkeLib </font>
︡bbc6b430-b6cf-4916-a85a-fb295b8b0951︡{"md":"<font size=5> Load ClarkeLib </font>"}︡
︠c6c6db02-1443-4b81-b799-f6a2147ec4cfs︠
load("clarkelib.sage")
︡730cd72e-faa6-4661-9869-aa2164dd275f︡
︠525cc65a-e025-4013-a9d8-e25bb53b1c66is︠
%md
<font size=5> Create a dictionary with tether material parameters </font>
︡59e90fa4-eb64-4253-9fc2-545694eb17dd︡{"md":"<font size=5> Create a dictionary with tether material parameters </font>"}︡
︠aeb18aa5-680f-455e-bc91-4c699971918es︠
tether_material = {'rho':1300, 'ss':50, 'ksafe':1.3}
︡8d4211b1-1455-4e3f-a61f-4007bdd54c13︡
︠62f65695-7d90-4d95-912c-f1637ae2774bis︠
%md
<font size=5> Create a Space Elevator object </font>
︡6f391ee8-8288-4773-9951-02a8cb1a8189︡{"md":"<font size=5> Create a Space Elevator object </font>"}︡
︠7b859e3e-8132-47df-b928-df8b1745d503s︠
EarthElevator = SpaceElevator(n_climbers = 7,
                              m_climber = 5e3,
                              anchor_force = 10e3,
                              material = tether_material,
                              planet = "Earth",
                              taper_func=('cubic', 'exp'),
                              khcw = 1.8,
                              anchor_safe='light')
︡8cac54f6-f4f0-41cd-8d09-46ebfcdb9631︡
︠011bc5c6-60a0-4188-84ba-bd3f65a81a21is︠
%md
<font size=5> Check accuracy of the equilibrium equation, N </font>
︡af650eb1-50ff-4d2e-b982-74aa180b3564︡{"md":"<font size=5> Check accuracy of the equilibrium equation, N </font>"}︡
︠f620cee8-844f-4231-a2d1-83b6208e038as︠
EarthElevator.tether.getEquilibriumAccuracy()
︡f7b9b055-eca0-4e4e-9486-462e30f1d5f1︡{"stdout":"1.74622982740402e-10\n"}︡
︠b00d5baf-5e3d-4cc3-9000-ef718d9cb9dais︠
%md
<font size=5> Display space elevator parameters </font>
︡168de7e1-359f-4eb4-8db4-9edd8252621c︡{"md":"<font size=5> Display space elevator parameters </font>"}︡
︠7eff0308-26da-4472-b727-80a4fa2e5916s︠
EarthElevator.plots.showData(width=10, height=7.1, fontsize=11)
︡2c4fc418-79ff-4b07-bad8-07887e175994︡{"once":false,"file":{"show":true,"uuid":"dbbe8c68-241c-45ac-81d7-4942775d0157","filename":"62addf41-8180-4076-bdfd-369f4c22ad64.svg"}}︡
︠7c2d7df5-b3de-430d-85c5-a1abd2432582is︠
%md
<font size=5> Analize space elevator mass dependency on the counterweight altitude </font>
︡256903f4-e5a2-4760-8da5-3e2990edce20︡{"md":"<font size=5> Analize space elevator mass dependency on the counterweight altitude </font>"}︡
︠65f1500e-f9c3-4a30-8dc4-ad53a04f7157s︠
hcw_start = EarthElevator.planet.HSYN*1.05                      # Initial counterweight altitude
hcw_end = EarthElevator.planet.HSYN*1.8                         # Final counterweight altitude
m0_min = EarthElevator.varyCWAltitude(hcw_start, hcw_end, 80)   # Run analysis
m0_min*1e-6                                                     # Print altitude of the lightest space elevator in the specified range of altitudes, Mm
︡0b29b8c2-9a64-49f7-b5a3-5e515eb58629︡{"stdout":"75.49937199502925\n"}︡
︠b8559a68-308b-4769-9bc9-460d061475d5is︠
%md
<font size=5> Plot analysis results </font>
︡462a6267-0ae0-4b3c-96ca-dadbb5b9375c︡{"md":"<font size=5> Plot analysis results </font>"}︡
︠aa1c090d-5e5a-4487-8f80-2f21fcbdf5d3s︠
EarthElevator.plots.mplplotMvsHcw(figsize = (6,8))
︡6d804324-1e00-42c1-b97d-e20ca4ac7ce0︡{"once":false,"file":{"show":true,"uuid":"196314de-8369-467a-a0b2-b1f278b0a234","filename":"dba83679-b0af-46c4-b6a1-499d0f8ba65e.svg"}}︡
︠b85b907c-0d9e-4203-9986-b2ec1947db17is︠
%md
<font size=5> Change counterweight altitude </font>
︡b9b7ae7a-5c24-45c2-9867-c2c21058a856︡{"md":"<font size=5> Change counterweight altitude </font>"}︡
︠5007ba4d-0316-4cd2-82ec-2ebca9252377s︠
EarthElevator.setHcw(100e6)
EarthElevator.plots.showData(width=10, height=7.1, fontsize=11)
︡ae8ee3d5-0e41-4cc7-987c-5e36a526e9b1︡{"once":false,"file":{"show":true,"uuid":"79c08aea-ab1a-425c-b6b7-c4d11304d5d5","filename":"7e086563-743a-4296-bba4-f3a0e7840ca1.svg"}}︡
︠74fa380b-d5d3-4f68-9a8e-4ac1bd1b52eais︠
%md
<font size=5> Analize tension in tether without climbers </font>
︡8f1e3abd-38ee-42b4-be74-b60746a28c5f︡{"md":"<font size=5> Analize tension in tether without climbers </font>"}︡
︠d094adef-8367-4ca4-af7f-f655b42defbbs︠
EarthElevator.plots.mplplotTetherProfile(climbers=False, figsize = (7,8))
︡128ac159-9c36-4975-909e-6b264f78bcc6︡{"once":false,"file":{"show":true,"uuid":"e44c12e4-8de2-4246-ae7c-20b578c9a409","filename":"501917cb-87cd-4eb2-89d7-d3029c7fce0d.svg"}}︡
︠cf371895-29cf-4788-ad2c-952fbc6b469cis︠
%md
<font size=5> Analize tension in tether with climbers </font>
︡cea588ab-02ca-4470-97c8-9ca8bf028cbf︡{"md":"<font size=5> Analize tension in tether with climbers </font>"}︡
︠f43c1b35-c0ab-4b29-9f17-0c07e1d2a3c7s︠
EarthElevator.plots.mplplotTetherProfile(climbers = True, figsize = (7,8))
︡2beaf73a-b4ef-4258-9843-d3524585f82f︡{"once":false,"file":{"show":true,"uuid":"317a8095-ee99-445c-9e98-11a899400e9c","filename":"62aae308-923c-4b04-85ad-957ff1f31998.svg"}}︡
︠358a04a1-0ae1-408e-9d07-107f86c1acf2is︠
%md
<font size=5> Now do the same for Mars Elevator </font>
︡959dc36a-3754-460c-b301-584fbd42d080︡{"md":"<font size=5> Now do the same for Mars Elevator </font>"}︡
︠6792ea1a-d343-4e4d-ba08-39c15423a0dbs︠
tether_material = {'rho':1300, 'ss':30, 'ksafe':1.3}
H_deimos = 23459000 - 3396*1e3
MarsElevator = SpaceElevator(n_climbers = 4,
                              m_climber = 5e3,
                              anchor_force = 10e3,
                              material = tether_material,
                              planet = "Mars",
                              taper_func=('cubic', 'exp'),
                              hcw = H_deimos*0.99,
                              anchor_safe='light')
H_deimos*1e-6
MarsElevator.plots.showData(width=10, height=7.1, fontsize=11)
︡0ae62373-c8b3-4c2b-a2d2-8f62c2aa4b32︡{"stdout":"20.0630000000000\n"}︡{"once":false,"file":{"show":true,"uuid":"c1de9087-0080-43c4-beaf-8e051addf5dd","filename":"4363b56f-17ff-4fe7-80ce-35d5945a6752.svg"}}︡
︠3ee78966-e14f-4598-8d22-c5ed1d916813︠









