from flask import Flask
from flask_bootstrap import Bootstrap
import flask_excel as excel

app = Flask(__name__,  static_url_path='', static_folder='static')

bootstrap = Bootstrap()
bootstrap.init_app(app)


excel.init_excel(app)

from app import routes
