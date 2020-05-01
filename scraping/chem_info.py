import urllib.request
from selenium import webdriver
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import time


def getCAS(MoleculeName, showImage=True):
    driver = webdriver.Chrome()
    # driver.minimize_window()
    driver.get("https://www.commonchemistry.org/")

    # input name
    driver.find_element_by_id("ctl07_TextSearchTerms").send_keys(MoleculeName)
    # search
    driver.find_element_by_id('ctl07_searchbutton').click()

    # select first item
    time.sleep(2)
    search_item_path = '//*[@id="searchresults"]/h2[1]/a'
    if len(driver.find_elements_by_xpath(search_item_path)) > 0:
        driver.find_element_by_xpath(search_item_path).click()

    # get cas number
    time.sleep(2)
    cas_path = '//*[@id="registryNumberLabel"]'
    if len(driver.find_elements_by_xpath(cas_path)) > 0:
        cas = driver.find_element_by_xpath(cas_path).text

        # image
        if showImage:
            image = driver.find_element_by_xpath('//*[@id="structureDiagram"]')
            image_url = image.get_attribute('src')
            png_name = '.\\temp\\'
            png_name += MoleculeName + '.png'
            urllib.request.urlretrieve(image_url, png_name)

            img = mpimg.imread(png_name)
            plt.imshow(img)
            plt.show()

        driver.close()
        return cas
    else:
        print('input molecule not existing.')
        driver.close()
        return ''


# getCAS('toluene')
end=1
