from manim import *

def dist(x1,y1,z1,x2,y2,z2):
    return ((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)**0.5

class Lattice:
    def __init__(self,a,b,c,numX,numY,numZ,x0,y0,z0, scene):
        self.a=a
        self.b=b
        self.c=c
        self.numX = numX
        self.numY = numY
        self.numZ = numZ
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.scene = scene

    def coordToPos(self,x,y,z):
        numX = self.numX - 1
        numY = self.numY - 1
        numZ = self.numZ - 1

        x_new = self.a*(-numX/2+x)+self.x0
        y_new = self.b*(-numY/2+y)+self.y0
        z_new = self.c*(-numZ/2+z)+self.z0
        return [x_new, y_new, z_new]

    def transform(self,x,y,z):
        #!!!
        self.scene.play(self.lattice_group.animate.move_to([x,y,z]))
        #self.lattice_group.move_to([x,y,z])
        self.x0 += x
        self.y0 += y
        self.z0 += z

    def drawLattice(self):
        dots = []
        

        x=0
        while x<self.numX:
            y=0
            while y<self.numY:
                z=0
                while z<self.numZ:

                    dots.append(Dot3D(point=self.coordToPos(x,y,z), radius=0.1))

                    z+=1
                y+=1
            x+=1

        lines = []
        y=0
        while y<self.numY:
            z=0
            while z<self.numZ:
                #line = Line(self.coordToPos(0, y, z), self.coordToPos(self.numX-1, y, z)).set_color(ORANGE)
                #lines.append(line)

                origin = self.coordToPos(0,0,self.numZ-1)
                line_length = self.numX - 1
                num_segments = 100
                maxDistance = dist(*self.coordToPos(self.numX-1,self.numY-1,self.numZ-1),*origin)
                colors = [PINK, BLACK]

                for i in range(num_segments):
                    start = np.array([line_length * i/num_segments, y, z])
                    end = np.array([line_length * (i+1)/num_segments, y, z])
                    midpoint = self.coordToPos(*((start + end) / 2))
                    distance = dist(*midpoint, *origin)
                    color = interpolate_color(colors[0], colors[-1], distance / maxDistance)


                    # start = np.array([i * line_length / num_segments, y, z])
                    # end = np.array([(i + 1) * line_length / num_segments, y, z])
                    # midpoint = (start + end) / 2
                    # distance = np.linalg.norm(midpoint - origin)
                    #color = interpolate_color(colors[0], colors[-1], distance / line_length)
                    lines.append(Line(self.coordToPos(*start), self.coordToPos(*end), color=color))

                z+=1
            y+=1
        lines = VGroup(*lines)
        
        p1 = self.coordToPos(0,0,0)
        p2 = self.coordToPos(1,0,0)
        p3 = self.coordToPos(0,1,0)
        p4 = self.coordToPos(0,0,1)

        # p2 = self.coordToPos(1,0,0)
        # p3 = self.coordToPos(0,1,0)
        # p4 = self.coordToPos(0,0,1)

        # b1 = BraceBetweenPoints(p2, p1)
        # b1text = Text('a').next_to(b1, DOWN)

        # b2 = BraceBetweenPoints(p1, p3)
        # b2text = Text('b').next_to(b2, RIGHT)

        # b3 = BraceBetweenPoints(p1, p4)
        # b3text = Text('c').next_to(b3, RIGHT)
        xVector = Arrow(p1, p2, buff=0, color=BLUE)
        yVector = Arrow(p1, p3, buff=0, color=GREEN)
        zVector = Arrow(p1, p4, buff=0, color=RED)

        dots = VGroup(*dots)

        # lattice_group = VGroup(lines,b1,b2,b1text,b2text,dots)
        lattice_group = VGroup(lines,dots,xVector,yVector,zVector)

        # dots.move_to(ORIGIN)
        #ThreeDScene.add(dots, b1, b2, b1text, b2text)
        #return lattice_group
        self.lattice_group = lattice_group
        self.scene.add(lattice_group)

    # def drawWave(self,x0,x1,y0,y1,z0,z1,n,A,t):
        

    #     x0,y0,z0=self.coordToPos(x0,y0,z0)
    #     x1,y1,z1=self.coordToPos(x1,y1,z1)

    #     # Calculate the distance and angle between the start and end points
    #     start_point = np.array([x0, y0, z0])
    #     end_point = np.array([x1, y1, z1])
    #     direction = end_point - start_point
    #     length = np.linalg.norm(direction)
    #     angle = np.arctan2(direction[1], direction[0])

    #     k = n*(2*PI/length)

    #     # Define the sine wave function along this direction
    #     def wave_func(t):
    #         x = t * length  # Parametrize x along the length of the line
    #         y = A * np.sin(k * x)  # Adjust amplitude and wave vector
    #         return start_point + np.array([x * np.cos(angle) - y * np.sin(angle),
    #                                        x * np.sin(angle) + y * np.cos(angle),
    #                                        0])

    #     # Create the sine wave as a parametric function
    #     wave = ParametricFunction(wave_func, t_range=[0, 1], color=BLUE)
        
    #     # Add the wave to the scene
    #     self.scene.add(wave)
    #     self.scene.play(Create(wave, run_time=t, rate_func=linear))

    def drawWave(self, x0, x1, y0, y1, z0, z1, n, A, t):
        # Convert coordinates to positions
        x0, y0, z0 = self.coordToPos(x0, y0, z0)
        x1, y1, z1 = self.coordToPos(x1, y1, z1)

        # Calculate the start and end points
        start_point = np.array([x0, y0, z0])
        end_point = np.array([x1, y1, z1])

        # Calculate the direction vector and normalize it
        direction = end_point - start_point
        length = np.linalg.norm(direction)
        unit_direction = direction / length

        # Calculate the wave number
        k = n * (2 * PI / length)

        # Create a vector orthogonal to the direction for the oscillation
        if np.allclose(unit_direction, [1, 0, 0]):
            orthogonal = np.array([0, 1, 0])
        else:
            orthogonal = np.cross(unit_direction, [1, 0, 0])
        orthogonal = orthogonal / np.linalg.norm(orthogonal)

        # Define the sine wave function along the direction vector
        def wave_func(t):
            x = t * length  # Parametrize x along the length of the line
            displacement = A * np.sin(k * x)
            return start_point + x * unit_direction + displacement * orthogonal

        # Create the sine wave as a 3D parametric function
        wave = ParametricFunction(wave_func, t_range=[0, 1], color=BLUE)

        # Add the wave to the scene
        self.scene.add(wave)
        self.scene.play(Create(wave, run_time=t, rate_func=linear))

        

class ReciprocalLattice(ThreeDScene):

    global textBoard 
    textBoard = {}

    def writeText(self,texts,positions,colors):
        h = lambda a: 3-(a-1)*1.25
        newTexts = []
        for text,pos,color in zip(texts,positions,colors):
            newText = MathTex(text,color=color).set_opacity(0)
            newText.move_to([3.75,h(pos),0])
            self.add_fixed_in_frame_mobjects(newText)
            textBoard[pos] = newText
            newTexts.append(newText)
        self.play(
            [text.animate.set_opacity(1) for text in newTexts]
        )

    def deleteText(self,positions):
        self.play(
            [textBoard[pos].animate.set_opacity(0) for pos in positions],
            run_time=0.75
        )
        for pos in positions:
            self.remove(textBoard[pos])
            textBoard[pos] = ""
    
    def construct(self):
        #self.next_section(skip_animations=True)

        lattice = Lattice(3,1.5,1,  3,3,3,  0,0,0, self)
        lattice.drawLattice()
        self.camera.set_focal_distance(1000)
        #self.move_camera(phi=60 * DEGREES, theta=-45 * DEGREES)
        self.camera.set_phi(60*DEGREES)
        self.camera.set_theta(-135*DEGREES)

        lattice.transform(-2.5,2,0.5)

        #  k = number of oscillations * (2pi/length)
        # lattice.drawWave(0.5,3,  1,2,  3,0.3,  0.5)
        #lattice.drawWave(1,2,1.5,  3,3,  0.3,  0.5)
        # lattice.drawWave(0,2.5,  3,-0.5,  7,0.3 ,  1)
        ##"\\rho(\\vec{r}) = \sum_{\\vec{R}} \delta(\\vec{r} - \\vec{R})",
        self.writeText(
            ["\\vec{R} = m \\vec{a}_1 + n \\vec{a}_2 + o \\vec{a}_3"],
            ## CHANGE COLORS OF LETTERS LATER ON
            [1],
            [YELLOW]
        )

        self.wait(20)

        self.deleteText([1])

        self.writeText(
            ["\Psi_k(\\vec{r}) = \Psi_0 \cdot e^{i \\vec{k} \cdot \\vec{r}}",
             "\Psi_k(\\vec{r})=\Psi_k(\\vec{r}+\\vec{R})",
             "\Psi_0 \cdot e^{i \\vec{k} \cdot \\vec{r}} = \Psi_0 \cdot e^{i \\vec{k} \cdot (\\vec{r} + \\vec{R})}",
            ],
            [1,2,3],
            [WHITE]*3
        )

        lattice.drawWave(0.5,1,  3,2,  3,2,   3,0.3,0.5)
        dot1 = Dot3D(point=lattice.coordToPos(1,2,2), radius=0.15, color=YELLOW)
        self.add(dot1)

        lattice.drawWave(1.5,2,  1,0,  1,0,   3,0.3,0.5)
        dot2 = Dot3D(point=lattice.coordToPos(2,0,0), radius=0.15, color=YELLOW)
        self.add(dot2)

        

        p1 = lattice.coordToPos(1,2,2)
        p2 = lattice.coordToPos(1,1,2)
        p3 = lattice.coordToPos(1,0,2)
        p4 = lattice.coordToPos(1,0,1)
        p5 = lattice.coordToPos(1,0,0)
        p6 = lattice.coordToPos(2,0,0)

        v1 = Arrow(p1, p2, buff=0, color=ORANGE)
        v2 = Arrow(p2, p3, buff=0, color=ORANGE)
        v3 = Arrow(p3, p4, buff=0, color=ORANGE)
        v4 = Arrow(p4, p5, buff=0, color=ORANGE)
        v5 = Arrow(p5, p6, buff=0, color=ORANGE)
        rgroup = VGroup(v1,v2,v3,v4,v5)

        self.play(FadeIn(v1, run_time=0.25))
        self.play(FadeIn(v2, run_time=0.25))
        self.play(FadeIn(v3, run_time=0.25))
        self.play(FadeIn(v4, run_time=0.25))
        self.play(FadeIn(v5, run_time=0.25))

        #self.next_section()

        v6 = Arrow(p1,p6,buff=0,color=ORANGE)
        self.play(
            ReplacementTransform(rgroup,v6)
        )

        #self.add(v1,v2,v3,v4,v5)

        ##ANIMATE WAVE AND VECTORS HERE
        self.writeText(
            ["e^{i\\vec{k} \cdot \\vec{R}} = 1",
             "\\vec{k} \cdot \\vec{R} = 2\pi l, \quad l \in \mathbb{Z}",
            ],
            [4,5],
            [WHITE,YELLOW]
        )
