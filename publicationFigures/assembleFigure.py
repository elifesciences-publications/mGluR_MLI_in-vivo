import svg_stack as ss
import os

doc = ss.Document()

layout1 = ss.HBoxLayout()
layout1.addSVG('WalkingFluorescenceBeforeAfterDrug_2.svg',alignment=ss.AlignTop|ss.AlignHCenter)
layout1.addSVGNoLayout('wheel_mouse_drug_application.svg',x=-1970,y=-40)
#layout1.addSVGNoLayout('spatial_illustration.svg',x=-1450,y=-35)

#layout2 = ss.VBoxLayout()

#layout2.addSVG('red_ball.svg',alignment=ss.AlignCenter)
#layout2.addSVG('red_ball.svg',alignment=ss.AlignCenter)
#layout2.addSVG('red_ball.svg',alignment=ss.AlignCenter)
#layout1.addLayout(layout2)

doc.setLayout(layout1)

figname = 'WalkingFluorescenceBeforeAfterDrugAssembled2'

doc.save(figname+'.svg')
os.system('inkscape -f '+str(figname)+'.svg -A '+str(figname)+'.pdf')
