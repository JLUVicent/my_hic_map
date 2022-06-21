class a():
    def __init__(self):
        self.b=0
        print(self.b)
        self.chengji()
        self.b=10
    def chengji(self):

        print(self.b*2)
if __name__=='__main__':
    value=a()
    print(value)