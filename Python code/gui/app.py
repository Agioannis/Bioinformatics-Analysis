from .tabs import EnhancedDNAToolGUI
import ttkbootstrap as tb

def main():
    root = tb.Window(themename='superhero')
    app = EnhancedDNAToolGUI(root)
    root.mainloop()

if __name__ == '__main__':
    main()
