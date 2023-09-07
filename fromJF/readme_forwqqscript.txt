
#####################################################
##############    ESTADO ACTUAL    ##################
#####################################################

Comentarios sobre el script selection_fromJF_wqq.py

Actualizar cada vez que se modifique el script. Ahora mismo se guardan, por comodidad, los root files para histogramas siempre en el mismo folder hist.../fromJF/wqq/ o ../fromJF/test/ pero esto 
puede llevar a confusion o a errores cuando se cambia el script, es por ello que en este documento se dejan los cambios realizados

-------------- 14/03/2023 ----------------

Ahora mismo se piden estos filtros:
     - Exactamente un muon bueno o electron bueno 
     - Al menos cuatro jets buenos, dos que no sean de los bottom
     - Dos btags, uno medium y otro loose
     - Corte en masa transversa de 50 GeV
     - Corte en masa invariante del dijet entre 30 y 330 GeV

Se ha comentado un filtro de cutbased para electrones que se estaba usando para las ntuplas en las que los electrones estaban sin ID

SF y otros pesos aplicados:
     - Están comentados los pesos que incluyen los SF de PUjetID ya que en las ntuplas usadas ahora no se pide el ID a los jets
     - El SF de ID para lepton se está poniendo y quitando para las pruebas en las que estamos ahora

Dónde se guardan los root files:
     - Normalmente se guardan en /nfs/cms/vazqueze/hists_ttbar/hists/higgs/fromJF/wqq/ pero ahora mismo en el script se guarda en ..../fromJF/test/
     - Cambiamos ahora a .../fromJF/test2/ para hacer dos pruebas simultáneas

-------------- 15/03/2023 ----------------

Ahora estamos corriendo cutbased, así que tomamos los datasets sin ID para el electrón y añadimos en los filtros el corte cutbased

SF y otros pesos aplicados:
     - Están comentados los pesos que incluyen los SF de PUjetID ya que en las ntuplas usadas ahora no se pide el ID a los jets
     - El SF de ID para lepton ahora mismo es el de cutbased
     - El SF de trig se quitó y se puso para varias pruebas pero ahora mismo está puesto

Dónde se guardan los root files:
     - En una carpeta .../fromJF/wqq/cutbased/

